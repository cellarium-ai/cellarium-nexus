"""Utility functions for submitting SOMA ingest and validation pipelines."""

import datetime
import os
import secrets

from django.conf import settings

from cellarium.nexus.backend.cell_management.models import OmicsDataset
from cellarium.nexus.backend.ingest_management.models import IngestSchema
from cellarium.nexus.backend.ingest_management.utils.soma_schema_utils import django_ingest_schema_to_pydantic
from cellarium.nexus.shared import schemas, utils
from cellarium.nexus.shared.utils import WorkspaceFileManager


def _generate_ingest_plan_path(*, base_name: str = "soma-ingests") -> str:
    """
    Generate timestamped ingest plan path for SOMA operations.

    :param base_name: Base directory name

    :return: GCS path for ingest plan
    """
    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    random_hex = secrets.token_hex(8)
    base_dir = (
        f"gs://{settings.BUCKET_NAME_PRIVATE}/{settings.BACKEND_PIPELINE_DIR}/{base_name}/{timestamp}_{random_hex}"
    )
    return f"{base_dir}/ingest_plan.json"


def _batch_uris(*, uris: list[str], batch_size: int) -> list[list[str]]:
    """
    Split URIs into batches of specified size.

    :param uris: List of URIs to batch
    :param batch_size: Number of URIs per batch

    :return: List of URI batches
    """
    batches = []
    for i in range(0, len(uris), batch_size):
        batches.append(uris[i : i + batch_size])
    return batches


def _generate_output_uris(*, input_uris: list[str], output_dir: str) -> list[str]:
    """
    Generate output URIs for sanitized files based on input filenames.

    :param input_uris: List of input GCS URIs
    :param output_dir: Base output directory

    :return: List of output URIs with same filenames in output_dir
    """
    output_dir = output_dir.rstrip("/")
    output_uris = []
    for uri in input_uris:
        filename = os.path.basename(uri)
        output_uris.append(f"{output_dir}/{filename}")
    return output_uris


def _dump_ingest_schema_to_gcs(*, ingest_schema: IngestSchema) -> str:
    """
    Dump Django IngestSchema to GCS as JSON and return the URI.

    :param ingest_schema: Django IngestSchema model instance

    :return: GCS URI to the dumped schema JSON file
    """
    pydantic_schema = django_ingest_schema_to_pydantic(django_schema=ingest_schema)

    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    schema_filename = f"ingest_schema_{timestamp}_{secrets.token_hex(6)}.json"

    bucket_name = settings.BUCKET_NAME_PRIVATE
    workspace_manager = WorkspaceFileManager(bucket_name=bucket_name)
    workspace_manager.save_json_to_bucket(
        data=pydantic_schema.model_dump(),
        remote_path=f"{settings.PIPELINE_CONFIGS_DIR}/{schema_filename}",
    )

    return f"gs://{bucket_name}/{settings.PIPELINE_CONFIGS_DIR}/{schema_filename}"


def create_soma_validate_sanitize_configs(
    *,
    input_h5ad_uris: list[str],
    ingest_schema_uri: str,
    output_base_dir: str,
    validation_report_id: int | None = None,
) -> tuple[list[schemas.component_configs.SomaValidateSanitizeConfig], list[str]]:
    """
    Create SomaValidateSanitizeConfig objects for batched validation.

    :param input_h5ad_uris: List of GCS paths to input h5ad files
    :param ingest_schema_uri: GCS URI to the dumped ingest schema JSON file
    :param output_base_dir: Base GCS directory for sanitized output files
    :param validation_report_id: Optional validation report ID for status tracking

    :return: Tuple of (list of configs, list of all output URIs)
    """
    all_output_uris = _generate_output_uris(input_uris=input_h5ad_uris, output_dir=output_base_dir)

    files_per_batch = settings.SOMA_FILES_PER_VALIDATION_BATCH
    batched_inputs = _batch_uris(uris=input_h5ad_uris, batch_size=files_per_batch)
    batched_outputs = _batch_uris(uris=all_output_uris, batch_size=files_per_batch)

    configs = []
    for input_batch, output_batch in zip(batched_inputs, batched_outputs):
        config = schemas.component_configs.SomaValidateSanitizeConfig(
            experiment_uri="",  # Not used by PreprocessingCoordinator for validation
            input_h5ad_uris=input_batch,
            output_h5ad_uris=output_batch,
            ingest_schema_uri=ingest_schema_uri,
            validation_report_id=validation_report_id,
            nexus_backend_api_url=settings.SITE_URL,
        )
        configs.append(config)

    return configs, all_output_uris


def create_soma_ingest_plan_config(
    *,
    sanitized_h5ad_uris: list[str],
    omics_dataset: OmicsDataset,
    ingest_schema_uri: str,
    ingest_plan_gcs_path: str,
) -> schemas.component_configs.SomaIngestPlanConfig:
    """
    Create SomaIngestPlanConfig for the ingest plan preparation step.

    :param sanitized_h5ad_uris: List of sanitized h5ad file URIs
    :param omics_dataset: OmicsDataset with schema and URI configured
    :param ingest_schema_uri: GCS URI to the dumped ingest schema JSON file
    :param ingest_plan_gcs_path: GCS path to store the ingest plan

    :raises ValueError: If omics_dataset has no schema or URI

    :return: SomaIngestPlanConfig
    """
    if not omics_dataset.schema:
        raise ValueError(f"OmicsDataset '{omics_dataset.name}' has no schema configured")
    if not omics_dataset.uri:
        raise ValueError(f"OmicsDataset '{omics_dataset.name}' has no URI configured")

    return schemas.component_configs.SomaIngestPlanConfig(
        experiment_uri=omics_dataset.uri,
        measurement_name="RNA",
        ingest_schema_uri=ingest_schema_uri,
        ingest_batch_size=settings.SOMA_FILES_PER_INGEST_PARTITION,
        h5ad_uris=sanitized_h5ad_uris,
        ingest_plan_gcs_path=ingest_plan_gcs_path,
        first_adata_gcs_path=sanitized_h5ad_uris[0] if sanitized_h5ad_uris else None,
    )


def create_soma_ingest_partition_configs(
    *,
    sanitized_h5ad_uris: list[str],
    omics_dataset: OmicsDataset,
    ingest_plan_gcs_path: str,
    parent_ingest_id: int,
) -> list[schemas.component_configs.SomaIngestPartitionConfig]:
    """
    Create SomaIngestPartitionConfig objects for each partition.

    :param sanitized_h5ad_uris: List of sanitized h5ad file URIs
    :param omics_dataset: OmicsDataset with URI configured
    :param ingest_plan_gcs_path: GCS path where ingest plan is stored
    :param parent_ingest_id: Parent Ingest record ID for tracking

    :return: List of partition configs, one per partition
    """
    if not omics_dataset.uri:
        raise ValueError(f"OmicsDataset '{omics_dataset.name}' has no URI configured")

    ingest_batch_size = settings.SOMA_FILES_PER_INGEST_PARTITION
    batched_uris = _batch_uris(uris=sanitized_h5ad_uris, batch_size=ingest_batch_size)
    num_partitions = len(batched_uris)

    configs = []
    for partition_index in range(num_partitions):
        config = schemas.component_configs.SomaIngestPartitionConfig(
            experiment_uri=omics_dataset.uri,
            ingest_plan_uri=ingest_plan_gcs_path,
            partition_index=partition_index,
            h5ad_file_paths=sanitized_h5ad_uris,
            nexus_backend_api_url=settings.SITE_URL,
            omics_dataset_name=omics_dataset.name,
            parent_ingest_id=parent_ingest_id,
        )
        configs.append(config)

    return configs


def submit_soma_validation_pipeline(
    *,
    input_h5ad_uris: list[str],
    ingest_schema: IngestSchema,
    output_directory_uri: str,
    validation_report_id: int,
) -> tuple[str, list[str]]:
    """
    Submit SOMA validation pipeline for h5ad files.

    Validate and sanitize files, output to specified directory.
    Create SomaValidateSanitizeConfig batches and submit
    run_soma_validate_sanitize_pipeline.

    :param input_h5ad_uris: List of GCS paths to input h5ad files
    :param ingest_schema: Django IngestSchema model for validation
    :param output_directory_uri: GCS directory for sanitized output files
    :param validation_report_id: Validation report ID for tracking

    :raises IOError: If config upload to GCS fails
    :raises google.api_core.exceptions.GoogleAPIError: If pipeline submission fails

    :return: Tuple of (pipeline_dashboard_url, list_of_output_sanitized_uris)
    """
    from cellarium.nexus.workflows.kubeflow.pipelines.soma import run_soma_validate_sanitize_pipeline
    from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline

    schema_uri = _dump_ingest_schema_to_gcs(ingest_schema=ingest_schema)

    configs, output_uris = create_soma_validate_sanitize_configs(
        input_h5ad_uris=input_h5ad_uris,
        ingest_schema_uri=schema_uri,
        output_base_dir=output_directory_uri,
        validation_report_id=validation_report_id,
    )

    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/{settings.PIPELINE_CONFIGS_DIR}"
    config_paths = utils.workflows_configs.dump_configs_to_bucket(configs=configs, bucket_path=configs_stage_dir)

    pipeline_url = submit_pipeline(
        pipeline_component=run_soma_validate_sanitize_pipeline,
        display_name=f"SOMA Validate/Sanitize - {ingest_schema.name}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={"validate_sanitize_configs": config_paths},
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
        labels={"application": settings.GCP_APPLICATION_BILLING_LABEL, "method": "soma_validation"},
    )

    return pipeline_url, output_uris


def submit_soma_ingest_pipeline(
    *,
    sanitized_h5ad_uris: list[str],
    omics_dataset: OmicsDataset,
    parent_ingest_id: int,
) -> str:
    """
    Submit SOMA ingest pipeline for pre-validated/sanitized h5ad files.

    Create SomaIngestPlanConfig and SomaIngestPartitionConfigs,
    then submit soma_ingest_pipeline.

    :param sanitized_h5ad_uris: List of GCS paths to sanitized h5ad files
    :param omics_dataset: OmicsDataset with schema and URI configured
    :param parent_ingest_id: Parent Ingest record ID for tracking

    :raises ValueError: If omics_dataset has no schema or URI, or URIs list is empty
    :raises IOError: If config upload to GCS fails
    :raises google.api_core.exceptions.GoogleAPIError: If pipeline submission fails

    :return: URL to the Vertex AI Pipeline dashboard
    """
    from cellarium.nexus.workflows.kubeflow.pipelines.soma import soma_ingest_pipeline
    from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline

    if not sanitized_h5ad_uris:
        raise ValueError("sanitized_h5ad_uris cannot be empty")

    ingest_plan_path = _generate_ingest_plan_path()
    schema_uri = _dump_ingest_schema_to_gcs(ingest_schema=omics_dataset.schema)

    ingest_plan_config = create_soma_ingest_plan_config(
        sanitized_h5ad_uris=sanitized_h5ad_uris,
        omics_dataset=omics_dataset,
        ingest_schema_uri=schema_uri,
        ingest_plan_gcs_path=ingest_plan_path,
    )

    partition_configs = create_soma_ingest_partition_configs(
        sanitized_h5ad_uris=sanitized_h5ad_uris,
        omics_dataset=omics_dataset,
        ingest_plan_gcs_path=ingest_plan_path,
        parent_ingest_id=parent_ingest_id,
    )

    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/{settings.PIPELINE_CONFIGS_DIR}"

    ingest_plan_config_paths = utils.workflows_configs.dump_configs_to_bucket(
        configs=[ingest_plan_config], bucket_path=configs_stage_dir
    )
    partition_config_paths = utils.workflows_configs.dump_configs_to_bucket(
        configs=partition_configs, bucket_path=configs_stage_dir
    )

    pipeline_url = submit_pipeline(
        pipeline_component=soma_ingest_pipeline,
        display_name=f"SOMA Ingest - {omics_dataset.name}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={
            "ingest_plan_config": ingest_plan_config_paths[0],
            "ingest_partition_configs": partition_config_paths,
        },
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
        labels={"application": settings.GCP_APPLICATION_BILLING_LABEL, "method": "soma_ingest"},
    )

    return pipeline_url
