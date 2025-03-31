from kfp import dsl
from cellarium.nexus.workflows.kubeflow import components


@dsl.pipeline(
    name="nexus-pipelines-ingest-data",
    description="Create ingest files and ingest data to BigQuery",
)
def ingest_data_pipeline(create_ingest_configs: list[str], ingest_config: str) -> None:
    """
    Create ingest files and then ingest them to BigQuery.
    
    :param create_ingest_configs: List of GCS paths to configuration files for create_ingest_files_job
    :param ingest_config: GCS path to configuration file for ingest_data_to_bigquery_job
    
    :raise: RuntimeError if any component fails
    """
    # Create a list to keep track of all create_ingest_files operations
    create_ingest_ops = []
    
    # First, execute all create_ingest_files operations in parallel
    with dsl.ParallelFor(items=create_ingest_configs, name="create-ingest-files-workers", parallelism=64) as item:
        create_ingest_op = components.create_ingest_files_job(gcs_config_path=item)
        create_ingest_ops.append(create_ingest_op)
    
    # Then, execute a single ingest_data_to_bigquery operation after all create_ingest_files operations are done
    ingest_op = components.ingest_data_to_bigquery_job(gcs_config_path=ingest_config)
    
    # Make sure all create_ingest_files operations are completed before starting the ingest operation
    for create_op in create_ingest_ops:
        ingest_op.after(create_op)


@dsl.pipeline(
    name="nexus-pipelines-extract-data",
    description="Extract data",
)
def extract_data_pipeline(prepare_extract_config: str, extract_configs: list[str]) -> None:
    """
    Extract data from BigQuery tables into AnnData files.
    
    :param prepare_extract_config: GCS path to configuration file for prepare_extract_job
    :param extract_configs: List of GCS paths to configuration files for extract_job
    
    :raise: RuntimeError if any component fails
    """
    # Run the prepare extract operation first
    prepare_extract_op = components.prepare_extract_job(gcs_config_path=prepare_extract_config)
    
    # Then run all extract operations in parallel
    with dsl.ParallelFor(items=extract_configs, name="extract-workers", parallelism=64) as item:
        extract_op = components.extract_job(gcs_config_path=item.extract_config)
        extract_op.after(prepare_extract_op)
