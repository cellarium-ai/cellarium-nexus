import math

from django.conf import settings
from google.cloud import bigquery

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.admin import constants
from cellarium.nexus.backend.cell_management.admin.utils import exceptions
from cellarium.nexus.omics_datastore.bq_ops import BigQueryDataOperator
from cellarium.nexus.shared import schemas
from cellarium.nexus.shared.utils import workflows_configs


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
    bq_data_operator = BigQueryDataOperator(
        client=bq_client, project=settings.GCP_PROJECT_ID, dataset=bigquery_dataset.name
    )
    return bq_data_operator.count_cells(filter_statements=filters)


def compose_extract_curriculum_configs(
    name: str,
    creator_id: int,
    feature_schema: models.FeatureSchema,
    bigquery_dataset: models.BigQueryDataset,
    extract_bin_size: int,
    categorical_column_count_limit: int,
    extract_bin_keys: list[str] | None = None,
    filters: dict | None = None,
    metadata_extra_columns: list[str] | None = None,
) -> tuple[schemas.component_configs.BQOpsPrepareExtract, list[schemas.component_configs.BQOpsExtract]]:
    """
    Compose extract curriculum configs for the Kubeflow pipeline.

    :param name: Prefix for extract table names
    :param creator_id: ID of a user who initiated the extract.
    :param feature_schema: Feature schema containing gene features
    :param bigquery_dataset: BigQuery dataset to extract from
    :param extract_bin_size: Number of cells per extract bin
    :param categorical_column_count_limit: Maximum number of categories per categorical column to be considered as
            categorical. If the number of categories exceeds this limit, the column will not be unified across all
            extract files.
    :param extract_bin_keys: Optional list of keys to use for binning the extract. If not provided, the keys will be
        assigned randomly.
    :param filters: Optional dictionary of filter statements to apply
    :param metadata_extra_columns: Optional list of additional metadata columns to include

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

    prepare_extract_config = schemas.component_configs.BQOpsPrepareExtract(
        name=name,
        project_id=settings.GCP_PROJECT_ID,
        nexus_backend_api_url=settings.SITE_URL,
        bigquery_dataset=bigquery_dataset.name,
        features=features,
        categorical_column_count_limit=categorical_column_count_limit,
        filters=filters,
        obs_columns=obs_columns,
        extract_bin_size=extract_bin_size,
        bucket_name=settings.BUCKET_NAME_PRIVATE,
        extract_bucket_path=extract_bucket_path,
        creator_id=creator_id,
        metadata_extra_columns=metadata_extra_columns,
        extract_bin_keys=extract_bin_keys,
    )

    extract_configs = []
    for start_bin in range(0, num_bins, constants.BINS_PER_WORKER):
        end_bin = min(start_bin + constants.BINS_PER_WORKER, num_bins)
        worker_bins = list(range(start_bin, end_bin))

        if not worker_bins:
            continue

        extract_obs_columns = obs_columns + metadata_extra_columns if metadata_extra_columns else obs_columns
        extract_configs.append(
            schemas.component_configs.BQOpsExtract(
                name=name,
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=bigquery_dataset.name,
                bins=worker_bins,
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                extract_bucket_path=extract_bucket_path,
                obs_columns=extract_obs_columns,
                max_workers=len(worker_bins),
            )
        )
    return prepare_extract_config, extract_configs


def compose_and_dump_configs(
    feature_schema: models.FeatureSchema,
    creator_id: int,
    bigquery_dataset: models.BigQueryDataset,
    name: str,
    extract_bin_size: int,
    categorical_column_count_limit: int,
    extract_bin_keys: list[str] | None = None,
    filters: dict | None = None,
    metadata_extra_columns: list[str] | None = None,
) -> tuple[str, list[str]]:
    """
    Compose extract pipeline configs and dump them to GCS bucket.

    Create prepare extract and extract configs based on the provided parameters,
    then save them to GCS and return their paths.

    :param feature_schema: Feature schema containing gene features to extract
    :param creator_id: ID of a user who initiated the extract.
    :param bigquery_dataset: BigQuery dataset to extract data from
    :param name: Name for the extract
    :param extract_bin_size: Number of cells per extract bin
    :param extract_bin_keys: Optional list of keys to use for binning the extract. If not provided, the keys will be
        assigned randomly.
    :param filters: Optional dictionary of filter statements to apply
    :param categorical_column_count_limit: Maximum number of categories per categorical column to be considered as
            categorical. If the number of categories exceeds this limit, the column will not be unified across all
            extract files.
    :param metadata_extra_columns: Optional list of additional metadata columns to include

    :raise ValueError: If extract_bin_size <= 0
    :raise exceptions.ZeroCellsReturnedError: If no cells match the filters
    :raise IOError: If there's an error writing configs to GCS
    :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

    :return: Tuple containing the prepare extract config path and a list of extract config paths
    """
    prepare_extract_config, extract_configs = compose_extract_curriculum_configs(
        feature_schema=feature_schema,
        creator_id=creator_id,
        bigquery_dataset=bigquery_dataset,
        name=name,
        extract_bin_size=extract_bin_size,
        categorical_column_count_limit=categorical_column_count_limit,
        extract_bin_keys=extract_bin_keys,
        filters=filters,
        metadata_extra_columns=metadata_extra_columns,
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
    creator_id: int,
    bigquery_dataset: models.BigQueryDataset,
    name: str,
    extract_bin_size: int,
    categorical_column_count_limit: int,
    extract_bin_keys: list[str] | None = None,
    filters: dict | None = None,
    metadata_extra_columns: list[str] | None = None,
) -> str:
    """
    Submit a Kubeflow pipeline for data extraction with the provided configurations.

    Generate extract configs, save them to GCS, and submit the pipeline for execution.

    :param feature_schema: Feature schema containing gene features to extract
    :param creator_id: ID of a user who initiated the extract.
    :param bigquery_dataset: BigQuery dataset to extract data from
    :param name: Name for the extract
    :param extract_bin_size: Number of cells per extract bin
    :param categorical_column_count_limit: Maximum number of categories per categorical column to be considered as
            categorical. If the number of categories exceeds this limit, the column will not be unified across all
            extract files.
    :param extract_bin_keys: Optional list of keys to use for binning the extract. If not provided, the keys will be
        assigned randomly.
    :param filters: Optional dictionary of filter statements to apply
    :param metadata_extra_columns: Optional list of additional metadata columns to include

    :raise ValueError: If extract_bin_size <= 0
    :raise exceptions.ZeroCellsReturnedError: If no cells match the filters
    :raise IOError: If there's an error writing configs to GCS
    :raise google.api_core.exceptions.GoogleAPIError: If pipeline submission fails

    :return: URL to the Vertex AI Pipeline dashboard for the submitted job
    """
    from cellarium.nexus.workflows.kubeflow.pipelines import extract_data_pipeline  # noqa
    from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline  # noqa

    prepare_extract_config_path, extract_config_paths = compose_and_dump_configs(
        feature_schema=feature_schema,
        creator_id=creator_id,
        bigquery_dataset=bigquery_dataset,
        name=name,
        categorical_column_count_limit=categorical_column_count_limit,
        extract_bin_size=extract_bin_size,
        extract_bin_keys=extract_bin_keys,
        filters=filters,
        metadata_extra_columns=metadata_extra_columns,
    )
    return submit_pipeline(
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
