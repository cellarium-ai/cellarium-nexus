import typing as t
from kfp import dsl
from cellarium.nexus.workflows.kubeflow import components


@dsl.pipeline(
    name="nexus-pipelines-create-ingest-files",
    description="Create ingest files in parallel",
)
def create_ingest_files_pipeline(create_ingest_configs: t.List[str]) -> None:
    """
    Create ingest files in parallel.

    :param create_ingest_configs: List of GCS paths to configuration files for create_ingest_files_job

    :raise: RuntimeError if any component fails
    """
    with dsl.ParallelFor(items=create_ingest_configs, name="create-ingest-files-workers", parallelism=64) as item:
        components.create_ingest_files_job(gcs_config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-ingest-to-bigquery",
    description="Ingest data to BigQuery",
)
def ingest_to_bigquery_pipeline(ingest_config: str) -> None:
    """
    Ingest data to BigQuery.

    :param ingest_config: GCS path to configuration file for ingest_data_to_bigquery_job

    :raise: RuntimeError if any component fails
    """
    components.ingest_data_to_bigquery_job(gcs_config_path=ingest_config)


@dsl.pipeline(
    name="nexus-pipelines-prepare-extract",
    description="Prepare tables for extraction",
)
def prepare_extract_pipeline(prepare_extract_config: str) -> None:
    """
    Prepare BigQuery tables for data extraction.

    :param prepare_extract_config: GCS path to configuration file for prepare_extract_job

    :raise: RuntimeError if any component fails
    """
    components.prepare_extract_job(gcs_config_path=prepare_extract_config)


@dsl.pipeline(
    name="nexus-pipelines-run-extracts",
    description="Run parallel extract operations",
)
def run_extracts_pipeline(extract_configs: t.List[str]) -> None:
    """
    Run multiple extract operations in parallel.

    :param extract_configs: List of GCS paths to configuration files for extract_job

    :raise: RuntimeError if any component fails
    """
    with dsl.ParallelFor(items=extract_configs, name="extract-workers", parallelism=64) as item:
        components.extract_job(gcs_config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-ingest-data",
    description="Create ingest files and ingest data to BigQuery",
)
def ingest_data_pipeline(create_ingest_configs: t.List[str], ingest_config: str) -> None:
    """
    Create ingest files and then ingest them to BigQuery.

    :param create_ingest_configs: List of GCS paths to configuration files for create_ingest_files_job
    :param ingest_config: GCS path to configuration file for ingest_data_to_bigquery_job

    :raise: RuntimeError if any component fails
    """
    # First run all create_ingest_files tasks in parallel
    create_op = create_ingest_files_pipeline(create_ingest_configs=create_ingest_configs)

    # Then run ingest task
    ingest_op = ingest_to_bigquery_pipeline(ingest_config=ingest_config)
    ingest_op.after(create_op)


@dsl.pipeline(
    name="nexus-pipelines-extract-data",
    description="Extract data",
)
def extract_data_pipeline(prepare_extract_config: str, extract_configs: t.List[str]) -> None:
    """
    Extract data from BigQuery tables into AnnData files.

    :param prepare_extract_config: GCS path to configuration file for prepare_extract_job
    :param extract_configs: List of GCS paths to configuration files for extract_job

    :raise: RuntimeError if any component fails
    """
    prepare_op = prepare_extract_pipeline(prepare_extract_config=prepare_extract_config)

    extract_op = run_extracts_pipeline(extract_configs=extract_configs)
    extract_op.after(prepare_op)
