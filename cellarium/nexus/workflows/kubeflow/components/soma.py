"""SOMA data operation components for Kubeflow pipelines."""

from cellarium.nexus.workflows.kubeflow import conf
from cellarium.nexus.workflows.kubeflow.utils import job

SOMA_EXTRACT_LABELS = {**conf.LABELS, "method": "soma_extract"}


@job.dsl_component_job(
    machine_type="e2-highmem-8",
    display_name="soma_extract",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_EXTRACT_LABELS,
)
def soma_extract_job(config_path: str):
    """
    Extract data from SOMA experiment into AnnData files.

    Read a SomaOpsExtract configuration from cloud storage and use
    SomaDataOpsCoordinator to extract data for the specified range indices.

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
        partition_index=params.partition_index,
        max_ranges_per_partition=params.max_ranges_per_partition,
        output_format=params.output_format,
        max_workers_extract=params.max_workers_extract,
        max_workers_shuffle=params.max_workers_shuffle,
    )


@job.dsl_component_job(
    machine_type="e2-standard-4",
    display_name="mark_soma_curriculum_as_finished",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=SOMA_EXTRACT_LABELS,
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
    )
