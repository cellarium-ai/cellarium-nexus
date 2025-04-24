import math

from django.conf import settings
from google.cloud import bigquery

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.admin import constants
from cellarium.nexus.backend.cell_management.admin.utils import exceptions
from cellarium.nexus.omics_datastore.bq_ops import bq_datastore_controller
from cellarium.nexus.shared import schemas
from cellarium.nexus.shared.utils import workflows_configs
from cellarium.nexus.workflows.kubeflow import component_configs
from cellarium.nexus.workflows.kubeflow.pipelines import extract_data_pipeline
from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline


def get_total_cell_in_bq_number(bigquery_dataset: models.BigQueryDataset, filters: dict | None = None) -> int:
    """
    Get total number of cells in BigQuery dataset that match the given filters.

    :param bigquery_dataset: BigQuery dataset containing cell data
    :param filters: Optional dictionary of filter statements to apply

    :raise: google.api_core.exceptions.GoogleAPIError: If BigQuery API request fails
    :raise: google.api_core.exceptions.ClientError: If client initialization fails

    :return: Total number of cells matching the filters
    """
    bq_client = bigquery.Client(project=settings.GCP_PROJECT_ID)
    controller = bq_datastore_controller.BQDatastoreController(
        client=bq_client, project=settings.GCP_PROJECT_ID, dataset=bigquery_dataset.name
    )
    return controller.count_cells(filter_statements=filters)


def compose_extract_curriculum_configs(
    name: str,
    feature_schema: models.FeatureSchema,
    bigquery_dataset: models.BigQueryDataset,
    extract_bin_size: int,
    creator_id: int,
    filters: dict | None = None,
) -> tuple[component_configs.BQOpsPrepareExtract, list[component_configs.BQOpsExtract]]:
    """
    Compose extract curriculum configs for the Kubeflow pipeline.

    :param name: Prefix for extract table names
    :param creator_id: ID of a user who initiated the extract.
    :param feature_schema: Feature schema containing gene features
    :param bigquery_dataset: BigQuery dataset to extract from
    :param extract_bin_size: Number of cells per extract bin
    :param filters: Optional dictionary of filter statements to apply

    :raise: ValueError: If extract_bin_size <= 0
    :raise: google.api_core.exceptions.GoogleAPIError: If BigQuery API request fails

    :return: Tuple of prepare extract config and list of extract configs
    """
    if extract_bin_size <= 0:
        raise ValueError(f"Extract bin size must be greater than 0. Received: {extract_bin_size}")

    total_cells = get_total_cell_in_bq_number(bigquery_dataset=bigquery_dataset, filters=filters)
    num_bins = math.ceil(total_cells / extract_bin_size)
    features = [
        schemas.FeatureSchema(id=idx, symbol=feature.symbol, ensemble_id=feature.ensemble_id)
        for idx, feature in enumerate(feature_schema.features.all())
    ]

    if total_cells == 0:
        raise exceptions.ZeroCellsReturnedError("BigQuery dataset contains no cells matching the filters.")

    obs_columns = constants.CELL_INFO_EXTRACT_COLUMNS

    extract_bucket_path = f"{settings.BACKEND_PIPELINE_DIR}/data-extracts/{name}"

    prepare_extract_config = component_configs.BQOpsPrepareExtract(
        name=name,
        project_id=settings.GCP_PROJECT_ID,
        nexus_backend_api_url=settings.SITE_URL,
        bigquery_dataset=bigquery_dataset.name,
        features=features,
        filters=filters,
        obs_columns=obs_columns,
        extract_bin_size=extract_bin_size,
        bucket_name=settings.BUCKET_NAME_PRIVATE,
        extract_bucket_path=extract_bucket_path,
        creator_id=creator_id,
    )

    extract_configs = []
    for start_bin in range(0, num_bins, constants.BINS_PER_WORKER):
        end_bin = min(start_bin + constants.BINS_PER_WORKER, num_bins)
        worker_bins = list(range(start_bin, end_bin))

        if not worker_bins:
            continue

        extract_configs.append(
            component_configs.BQOpsExtract(
                name=name,
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=bigquery_dataset.name,
                bins=worker_bins,
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                extract_bucket_path=extract_bucket_path,
                obs_columns=obs_columns,
                max_workers=len(worker_bins),
            )
        )
    return prepare_extract_config, extract_configs


def compose_and_dump_configs(
    feature_schema: models.FeatureSchema,
    bigquery_dataset: models.BigQueryDataset,
    extract_table_prefix: str,
    extract_bin_size: int,
    filters: dict | None = None,
) -> tuple[str, list[str]]:
    prepare_extract_config, extract_configs = compose_extract_curriculum_configs(
        feature_schema=feature_schema,
        bigquery_dataset=bigquery_dataset,
        name=extract_table_prefix,
        extract_bin_size=extract_bin_size,
        filters=filters,
    )
    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"

    prepare_extract_config_path = workflows_configs.dump_configs_to_bucket(
        configs=[prepare_extract_config], bucket_path=configs_stage_dir
    )[0]
    extract_config_paths = workflows_configs.dump_configs_to_bucket(
        configs=extract_configs, bucket_path=configs_stage_dir
    )

    return prepare_extract_config_path, extract_config_paths


def submit_extract_pipeline(
    feature_schema: models.FeatureSchema,
    bigquery_dataset: models.BigQueryDataset,
    extract_table_prefix: str,
    extract_bin_size: int,
    filters: dict | None = None,
) -> None:
    prepare_extract_config_path, extract_config_paths = compose_and_dump_configs(
        feature_schema=feature_schema,
        bigquery_dataset=bigquery_dataset,
        extract_table_prefix=extract_table_prefix,
        extract_bin_size=extract_bin_size,
        filters=filters,
    )
    submit_pipeline(
        pipeline_component=extract_data_pipeline,
        display_name=f"Nexus Extract Data - {bigquery_dataset.name}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={
            "prepare_extract_config": prepare_extract_config_path,
            "extract_configs": extract_config_paths,
        },
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
    )
