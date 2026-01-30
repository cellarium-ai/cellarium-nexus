"""SOMA data operation components for Kubeflow pipelines."""

from cellarium.nexus.workflows.kubeflow import conf
from cellarium.nexus.workflows.kubeflow.utils import job

SOMA_RANDOMIZED_EXTRACT_LABELS = {**conf.LABELS, "method": "soma_randomized_extract"}
SOMA_GROUPED_EXTRACT_LABELS = {**conf.LABELS, "method": "soma_grouped_extract"}
SOMA_INGEST_LABELS = {**conf.LABELS, "method": "soma_ingest"}


@job.dsl_component_job(
    machine_type="e2-standard-8",
    display_name="soma_randomized_extract",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_RANDOMIZED_EXTRACT_LABELS,
)
def soma_randomized_extract_job(config_path: str):
    """
    Extract data from SOMA experiment using randomized extraction with shuffle.

    Read a SomaOpsExtract configuration from cloud storage and use
    SomaDataOpsCoordinator to extract data with cell randomization.

    :param config_path: Path to the configuration file in cloud storage

    :raise Exception: If extraction fails
    """
    from cellarium.nexus.coordinator import SomaDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import SomaOpsExtract

    params = utils.workflows_configs.read_component_config(gcs_path=config_path, schema_class=SomaOpsExtract)

    coordinator = SomaDataOpsCoordinator(
        experiment_uri=params.experiment_uri,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bucket_name=params.bucket_name,
    )
    coordinator.run_soma_extract(
        extract_name=params.extract_name,
        curriculum_metadata_path=params.extract_metadata_path,
        extract_bucket_path=params.extract_bucket_path,
        extract_type="randomized",
        partition_index=params.partition_index,
        output_format=params.output_format,
        max_workers_extract=params.max_workers_extract,
        max_ranges_per_partition=params.max_ranges_per_partition,
        max_workers_shuffle=params.max_workers_shuffle,
    )


@job.dsl_component_job(
    machine_type="e2-standard-16",
    display_name="soma_grouped_extract",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_GROUPED_EXTRACT_LABELS,
)
def soma_grouped_extract_job(config_path: str):
    """
    Extract data from SOMA experiment using grouped extraction.

    Read a SomaOpsExtract configuration from cloud storage and use
    SomaDataOpsCoordinator to extract data keeping cells from the same group together.

    :param config_path: Path to the configuration file in cloud storage

    :raise Exception: If extraction fails
    """
    from cellarium.nexus.coordinator import SomaDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import SomaOpsExtract

    params = utils.workflows_configs.read_component_config(gcs_path=config_path, schema_class=SomaOpsExtract)

    coordinator = SomaDataOpsCoordinator(
        experiment_uri=params.experiment_uri,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bucket_name=params.bucket_name,
    )
    coordinator.run_soma_extract(
        extract_name=params.extract_name,
        curriculum_metadata_path=params.extract_metadata_path,
        extract_bucket_path=params.extract_bucket_path,
        extract_type="grouped",
        partition_index=params.partition_index,
        output_format=params.output_format,
        max_workers_extract=params.max_workers_extract,
        max_bins_per_partition=params.max_bins_per_partition,
    )


@job.dsl_component_job(
    machine_type="e2-standard-4",
    display_name="mark_soma_curriculum_as_finished",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_RANDOMIZED_EXTRACT_LABELS,
)
def mark_soma_curriculum_as_finished_job(config_path: str):
    """
    Mark a SOMA curriculum as finished and update with final metadata.

    Read a SomaOpsExtract configuration from cloud storage and use
    SomaDataOpsCoordinator to update the curriculum status to SUCCEEDED.

    :param config_path: Path to the configuration file in cloud storage

    :raise Exception: If any error occurs during the process
    """
    from cellarium.nexus.coordinator import SomaDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import SomaOpsExtract

    params = utils.workflows_configs.read_component_config(gcs_path=config_path, schema_class=SomaOpsExtract)

    coordinator = SomaDataOpsCoordinator(
        experiment_uri=params.experiment_uri,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bucket_name=params.bucket_name,
    )
    coordinator.mark_soma_curriculum_as_finished(
        extract_name=params.extract_name,
        curriculum_metadata_path=params.extract_metadata_path,
        extract_bucket_path=params.extract_bucket_path,
        extract_type=params.extract_type,
    )


@job.dsl_component_job(
    machine_type="e2-highmem-8",
    display_name="soma_validate_sanitize",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_INGEST_LABELS,
)
def soma_validate_sanitize_job(config_path: str):
    """
    Validate and sanitize h5ad files for SOMA ingest.

    Read a SomaValidateSanitizeConfig from GCS and use SomaIngestCoordinator
    to validate and sanitize files for ingestion.
    """
    from cellarium.nexus.coordinator import SomaIngestCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import SomaValidateSanitizeConfig

    params = utils.workflows_configs.read_component_config(
        gcs_path=config_path, schema_class=SomaValidateSanitizeConfig
    )

    coordinator = SomaIngestCoordinator(experiment_uri=params.experiment_uri)
    coordinator.validate_and_sanitize_files(
        input_h5ad_uris=params.input_h5ad_uris,
        output_h5ad_uris=params.output_h5ad_uris,
        ingest_schema=params.ingest_schema,
        max_bytes_per_file=params.max_bytes_per_file,
    )


@job.dsl_component_job(
    machine_type="e2-standard-4",
    display_name="soma_prepare_ingest_plan",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_INGEST_LABELS,
)
def soma_prepare_ingest_plan_job(config_path: str):
    """
    Prepare a SOMA ingest plan and write it to GCS.

    Read a SomaIngestPlanConfig from GCS and use SomaIngestCoordinator
    to compute and save the ingest plan.
    """
    from cellarium.nexus.coordinator import SomaIngestCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import SomaIngestPlanConfig

    params = utils.workflows_configs.read_component_config(gcs_path=config_path, schema_class=SomaIngestPlanConfig)

    coordinator = SomaIngestCoordinator(experiment_uri=params.experiment_uri)
    ingest_plan = coordinator.prepare_ingest_plan(
        h5ad_uris=params.h5ad_uris,
        measurement_name=params.measurement_name,
        ingest_schema=params.ingest_schema,
        ingest_batch_size=params.ingest_batch_size,
        first_adata_gcs_path=params.first_adata_gcs_path,
    )
    coordinator.save_ingest_plan_to_gcs(ingest_plan=ingest_plan, ingest_plan_gcs_path=params.ingest_plan_gcs_path)


@job.dsl_component_job(
    machine_type="e2-highmem-8",
    display_name="soma_ingest_partition",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_INGEST_LABELS,
)
def soma_ingest_partition_job(config_path: str):
    """
    Ingest a SOMA partition.

    Read a SomaIngestPartitionConfig from GCS and use SomaIngestCoordinator
    to ingest a single partition of h5ad files.
    """
    from cellarium.nexus.coordinator import SomaIngestCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import SomaIngestPartitionConfig

    params = utils.workflows_configs.read_component_config(gcs_path=config_path, schema_class=SomaIngestPartitionConfig)

    coordinator = SomaIngestCoordinator(experiment_uri=params.experiment_uri)
    ingest_plan = coordinator.load_ingest_plan_from_gcs(ingest_plan_gcs_path=params.ingest_plan_gcs_path)
    coordinator.ingest_partition(
        ingest_plan=ingest_plan,
        partition_index=params.partition_index,
    )
