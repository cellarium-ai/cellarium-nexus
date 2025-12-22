"""Validation components for Kubeflow pipelines."""

from cellarium.nexus.workflows.kubeflow import conf
from cellarium.nexus.workflows.kubeflow.utils import job


@job.dsl_component_job(
    machine_type="e2-highmem-4",
    display_name="validate_anndata_files",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=conf.LABELS,
)
def validate_anndata_files_job(gcs_config_path: str):
    """
    Validate multiple AnnData files and report validation results.

    Download each AnnData file from GCS, apply validation methods, and report
    results to the Nexus backend API.

    :param gcs_config_path: Path to the configuration file in GCS

    :raise ValidationError: If validation fails for any file
    """
    from cellarium.nexus.coordinator import NexusDataValidationCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import ValidationConfig

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=ValidationConfig)

    coordinator = NexusDataValidationCoordinator(
        nexus_backend_api_url=params.nexus_backend_api_url,
        validation_report_id=params.validation_report_id,
        max_bytes_valid_per_file=params.max_bytes_valid_per_file,
    )
    coordinator.validate_and_report_multiple(
        adata_gcs_paths=params.adata_gcs_paths,
        validation_methods=params.validation_methods,
    )
