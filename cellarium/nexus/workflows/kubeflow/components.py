from cellarium.nexus.workflows.kubeflow import conf
from cellarium.nexus.workflows.kubeflow.utils import job

INGEST_LABELS = {**conf.LABELS, "method": "ingest"}
EXTRACT_LABELS = {**conf.LABELS, "method": "extract"}


@job.dsl_component_job(
    machine_type="e2-highmem-8",
    display_name="create_ingest_files",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=INGEST_LABELS,
)
def create_ingest_files_job(gcs_config_path: str):
    """
    Create ingest files from input data sources.

    Read a CreateIngestFilesConfig from GCS and use BQDataOpsCoordinator to create
    ingest files for the specified data source.

    :param gcs_config_path: GCS path to the configuration file

    :raise Exception: If file creation fails
    """
    from cellarium.nexus.coordinator import BQDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import CreateIngestFilesConfig

    params = utils.workflows_configs.read_component_config(
        gcs_path=gcs_config_path, schema_class=CreateIngestFilesConfig
    )

    coordinator = BQDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.create_ingest_files(
        input_file_path=params.input_file_path,
        tag=params.tag,
        bigquery_dataset=params.bigquery_dataset,
        bucket_name=params.bucket_name,
        bucket_stage_dir=params.bucket_stage_dir,
        column_mapping=params.column_mapping,
        max_input_data_size=params.max_input_data_size,
        uns_keys_to_keep=params.uns_keys_to_keep,
        validation_methods=params.validation_methods,
    )


@job.dsl_component_job(
    machine_type="e2-standard-4",
    display_name="ingest_data_to_bigquery",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=INGEST_LABELS,
)
def ingest_data_to_bigquery_job(gcs_config_path: str):
    """
    Ingest prepared data into BigQuery tables.

    Read an IngestFilesConfig from GCS and use BQDataOpsCoordinator to ingest
    data from the specified bucket paths into BigQuery tables.

    :param gcs_config_path: GCS path to the configuration file

    :raise Exception: If ingestion fails
    """
    from cellarium.nexus.coordinator import BQDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import IngestFilesConfig

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=IngestFilesConfig)

    coordinator = BQDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.ingest_data_to_bigquery_parallel(
        bucket_name=params.bucket_name,
        bucket_stage_dirs=params.bucket_stage_dirs,
        num_workers=params.num_workers,
    )


@job.dsl_component_job(
    machine_type="e2-standard-4",
    display_name="prepare_extract",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=EXTRACT_LABELS,
)
def prepare_extract_job(gcs_config_path: str):
    """
    Prepare BigQuery tables for data extraction.

    Read a BQOpsPrepareExtract configuration from GCS and use BQDataOpsCoordinator
    to prepare extract tables based on the specified features and filters.

    :param gcs_config_path: GCS path to the configuration file

    :raise Exception: If preparation fails
    """
    from cellarium.nexus.coordinator import BQDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import BQOpsPrepareExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    coordinator = BQDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.prepare_extract_tables(
        extract_name=params.extract_name,
        features=params.features,
        categorical_column_count_limit=params.categorical_column_count_limit,
        creator_id=params.creator_id,
        filters=params.filters,
        extract_bin_size=params.extract_bin_size,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
        obs_columns=params.obs_columns,
        metadata_extra_columns=params.metadata_extra_columns,
        extract_bin_keys=params.extract_bin_keys,
    )


@job.dsl_component_job(
    machine_type="e2-standard-32",
    display_name="extract",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=EXTRACT_LABELS,
)
def extract_job(gcs_config_path: str):
    """
    Extract data from prepared tables into AnnData files.

    Read a BQOpsExtract configuration from GCS and use BQDataOpsCoordinator
    to extract data from BigQuery tables into AnnData files for the specified bins.

    :param gcs_config_path: GCS path to the configuration file

    :raise Exception: If extraction fails
    """
    from cellarium.nexus.coordinator import BQDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import BQOpsExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsExtract)

    coordinator = BQDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.extract_data(
        extract_name=params.extract_name,
        bins=params.bins,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
        obs_columns=params.obs_columns,
        max_workers=params.max_workers,
    )


@job.dsl_component_job(
    machine_type="e2-standard-4",
    display_name="nexus_mark_curriculum_as_finished",
    base_image=conf.BASE_IMAGE,
    service_account=conf.SERVICE_ACCOUNT,
    labels=EXTRACT_LABELS,
)
def mark_curriculum_as_finished_job(gcs_config_path: str):
    """
    Mark a curriculum as finished and succeeded, including extract metadata.

    Retrieves metadata information from the metadata file and updates the curriculum
    with this information, including the number of bins, extract files path, and
    metadata file path.

    :param gcs_config_path: Path to the configuration file in GCS

    :raise Exception: If any error occurs during the process
    """
    from cellarium.nexus.coordinator import BQDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.shared.schemas.component_configs import BQOpsPrepareExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    coordinator = BQDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.mark_curriculum_as_finished(
        extract_name=params.extract_name,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
    )


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

    Downloads each AnnData file from GCS, applies validation methods, and reports results to the Nexus backend API.

    :param gcs_config_path: Path to the configuration file in GCS

    :raise: ValidationError if validation fails for any file
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
