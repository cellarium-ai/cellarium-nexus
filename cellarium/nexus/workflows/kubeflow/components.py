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
    from cellarium.nexus.nexus_data_controller import NexusDataController
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import IngestTaskConfig

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=IngestTaskConfig)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.create_ingest_files(
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
    from cellarium.nexus.nexus_data_controller import NexusDataController
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import IngestTaskConfig

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=IngestTaskConfig)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.ingest_data_to_bigquery(
        bucket_name=params.bucket_name,
        bucket_stage_dir=params.ingest_bucket_path,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE, machine_type="e2-standard-4", display_name="prepare_extract", service_account=SERVICE_ACCOUNT
)
def prepare_extract_job(gcs_config_path: str):
    from cellarium.nexus.nexus_data_controller import NexusDataController
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsPrepareExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.prepare_extract_tables(
        extract_name=params.name,
        features=params.features,
        creator_id=params.creator_id,
        filters=params.filters,
        obs_columns=params.obs_columns,
        extract_bin_size=params.extract_bin_size,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
    )


@job.dsl_component_job(
    base_image=BASE_IMAGE, machine_type="e2-standard-32", display_name="extract", service_account=SERVICE_ACCOUNT
)
def extract_job(gcs_config_path: str):
    from cellarium.nexus.nexus_data_controller import NexusDataController
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsExtract)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.extract_data(
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
    display_name="mark_curriculum_as_finished",
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
    from cellarium.nexus.nexus_data_controller import NexusDataController
    from cellarium.nexus.shared import utils
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsPrepareExtract

    params = utils.workflows_configs.read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.mark_curriculum_as_finished(
        extract_name=params.name,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
    )
