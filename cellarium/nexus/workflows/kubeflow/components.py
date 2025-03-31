from cellarium.nexus.workflows.kubeflow.utils import job


@job.dsl_component_job(base_image="python:3.10", machine_type="e2-standard-4", display_name="create_ingest_files")
def create_ingest_files_job(gcs_config_path: str):
    from cellarium.nexus.workflows.kubeflow.utils.component_configs import read_component_config
    from cellarium.nexus.workflows.kubeflow.component_configs import CreateIngestFiles
    from cellarium.nexus.omics_datastore.controller import NexusDataController

    params = read_component_config(gcs_path=gcs_config_path, schema_class=CreateIngestFiles)

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


@job.dsl_component_job(base_image="python:3.10", machine_type="e2-standard-16", display_name="ingest_data_to_bigquery")
def ingest_data_to_bigquery_job(gcs_config_path: str):
    from cellarium.nexus.workflows.kubeflow.utils.component_configs import read_component_config
    from cellarium.nexus.workflows.kubeflow.component_configs import IngestDataToBigQuery
    from cellarium.nexus.omics_datastore.controller import NexusDataController

    params = read_component_config(gcs_path=gcs_config_path, schema_class=IngestDataToBigQuery)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.ingest_data_to_bigquery(
        bucket_name=params.bucket_name,
        bucket_stage_dir=params.ingest_bucket_path,
    )


@job.dsl_component_job(base_image="python:3.10", machine_type="e2-standard-4", display_name="prepare_extract")
def prepare_extract_job(gcs_config_path: str):
    from cellarium.nexus.workflows.kubeflow.utils.component_configs import read_component_config
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsPrepareExtract
    from cellarium.nexus.omics_datastore.controller import NexusDataController

    params = read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsPrepareExtract)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.prepare_extract_tables(
        extract_table_prefix=params.extract_table_prefix,
        features=params.features,
        filters=params.filters,
        obs_columns=params.obs_columns,
        extract_bin_size=params.extract_bin_size,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
    )


@job.dsl_component_job(base_image="python:3.10", machine_type="e2-standard-32", display_name="extract")
def extract_job(gcs_config_path: str):
    from cellarium.nexus.workflows.kubeflow.utils.component_configs import read_component_config
    from cellarium.nexus.workflows.kubeflow.component_configs import BQOpsExtract
    from cellarium.nexus.omics_datastore.controller import NexusDataController

    params = read_component_config(gcs_path=gcs_config_path, schema_class=BQOpsExtract)

    controller = NexusDataController(
        project_id=params.project_id,
        nexus_backend_api_url=params.nexus_backend_api_url,
        bigquery_dataset=params.bigquery_dataset,
    )
    controller.extract_data(
        extract_table_prefix=params.extract_table_prefix,
        bins=params.bins,
        bucket_name=params.bucket_name,
        extract_bucket_path=params.extract_bucket_path,
        obs_columns=params.obs_columns,
        max_workers=params.max_workers,
    )
