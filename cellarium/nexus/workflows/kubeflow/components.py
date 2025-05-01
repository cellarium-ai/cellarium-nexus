from cellarium.nexus.workflows.kubeflow.utils import job

BASE_IMAGE = "us-central1-docker.pkg.dev/dsp-cellarium/nexus/nexus-workflows:latest"
SERVICE_ACCOUNT = "vertex-pipelines-sa@dsp-cellarium.iam.gserviceaccount.com"


@job.dsl_component_job(
    base_image=BASE_IMAGE,
    machine_type="e2-standard-8",
    display_name="create_ingest_files",
    service_account=SERVICE_ACCOUNT,
)
def create_ingest_files_job(gcs_config_path: str):
    """
    Create ingest files from input data sources.

    Reads a CreateIngestFilesConfig from GCS and uses NexusDataOpsCoordinator to create
    ingest files for the specified data source.

    :param gcs_config_path: GCS path to the configuration file

    :raise: Exception: If file creation fails
    """
    from cellarium.nexus.coordinator import NexusDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import CreateIngestFilesConfig

    params = utils.workflows_configs.read_component_config(
        gcs_path=gcs_config_path, schema_class=CreateIngestFilesConfig
    )

    coordinator = NexusDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.create_ingest_files(
        input_file_path=params.data_source_path,
        tag=params.tag,
        bigquery_dataset=params.bigquery_dataset,
        bucket_name=params.bucket_name,
        bucket_stage_dir=params.ingest_bucket_path,
        column_mapping=params.metadata_columns,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE,
    machine_type="e2-standard-4",
    display_name="ingest_data_to_bigquery",
    service_account=SERVICE_ACCOUNT,
)
def ingest_data_to_bigquery_job(gcs_config_path: str):
    """
    Ingest prepared data into BigQuery tables.

    Reads an IngestFilesConfig from GCS and uses NexusDataOpsCoordinator to ingest
    data from the specified bucket paths into BigQuery tables.

    :param gcs_config_path: GCS path to the configuration file

    :raise: Exception: If ingestion fails
    """
    from cellarium.nexus.coordinator import NexusDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import IngestFilesConfig

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=IngestFilesConfig)

    coordinator = NexusDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.ingest_data_to_bigquery_parallel(
        bucket_name=params.bucket_name,
        bucket_stage_dirs=params.ingest_bucket_paths,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE, machine_type="e2-standard-4", display_name="prepare_extract", service_account=SERVICE_ACCOUNT
)
def prepare_extract_job(gcs_config_path: str):
    """
    Prepare BigQuery tables for data extraction.

    Reads a BQOpsPrepareExtract configuration from GCS and uses NexusDataOpsCoordinator
    to prepare extract tables based on the specified features and filters.

    :param gcs_config_path: GCS path to the configuration file

    :raise: Exception: If preparation fails
    """
    from cellarium.nexus.coordinator import NexusDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsPrepareExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    coordinator = NexusDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.prepare_extract_tables(
        extract_name=params.name,
        features=params.features,
        categorical_column_count_limit=params.categorical_column_count_limit,
        creator_id=params.creator_id,
        filters=params.filters,
        extract_bin_size=params.extract_bin_size,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
        obs_columns=params.obs_columns,
        metadata_extra_columns=params.metadata_extra_columns,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE, machine_type="e2-standard-32", display_name="extract", service_account=SERVICE_ACCOUNT
)
def extract_job(gcs_config_path: str):
    """
    Extract data from prepared tables into AnnData files.

    Reads a BQOpsExtract configuration from GCS and uses NexusDataOpsCoordinator
    to extract data from BigQuery tables into AnnData files for the specified bins.

    :param gcs_config_path: GCS path to the configuration file

    :raise: Exception: If extraction fails
    """
    from cellarium.nexus.coordinator import NexusDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsExtract)

    coordinator = NexusDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.extract_data(
        extract_name=params.name,
        bins=params.bins,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
        obs_columns=params.obs_columns,
        max_workers=params.max_workers,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE,
    machine_type="e2-standard-4",
    display_name="nexus_mark_curriculum_as_finished",
    service_account=SERVICE_ACCOUNT,
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
    from cellarium.nexus.coordinator import NexusDataOpsCoordinator
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsPrepareExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    coordinator = NexusDataOpsCoordinator(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    coordinator.mark_curriculum_as_finished(
        extract_name=params.name,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE,
    machine_type="e2-standard-4",
    display_name="validate_anndata_files",
    service_account=SERVICE_ACCOUNT,
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
    from cellarium.nexus.workflows.kubeflow.component_configs import ValidationConfig

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
