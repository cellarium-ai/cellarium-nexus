import typing as t

from kfp import dsl

from cellarium.nexus.workflows.kubeflow import components


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
    name="nexus-pipelines-create-ingest-data",
    description="For each input config, create ingest files and then ingest data to BigQuery, in parallel.",
)
def create_ingest_files_parallel_pipeline(create_ingest_files_configs: t.List[str]) -> None:
    """
    For each input task configuration, create ingest files and then ingest them to BigQuery.
    These create->ingest sequences run in parallel for all input configurations.

    :param create_ingest_files_configs: List of GCS paths to the combined IngestTaskConfig files.

    :raise: RuntimeError if any component fails
    """
    with dsl.ParallelFor(items=create_ingest_files_configs, name="create-ingest-files-workers", parallelism=64) as item:
        components.create_ingest_files_job(gcs_config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-ingest-data-to-bigquery",
    description="Ingest data to BigQuery",
)
def ingest_data_to_bigquery_pipeline(create_ingest_files_configs: t.List[str], ingest_task_config: str) -> None:
    create_ingest_files_op = create_ingest_files_parallel_pipeline(
        create_ingest_files_configs=create_ingest_files_configs
    )
    ingest_data_to_bigquery_op = components.ingest_data_to_bigquery_job(gcs_config_path=ingest_task_config)
    ingest_data_to_bigquery_op.after(create_ingest_files_op)


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
    mark_curriculum_as_finished_op = components.mark_curriculum_as_finished_job(gcs_config_path=prepare_extract_config)
    mark_curriculum_as_finished_op.after(extract_op)


@dsl.pipeline(
    name="nexus-pipelines-validate-ingest-files",
    description="Validate ingest AnnData files and report results",
)
def validate_anndata_pipeline(validation_configs: t.List[str]) -> None:
    """
    Validate multiple AnnData files and report validation results.

    :param validation_configs: GCS path to configuration file for validate_anndata_files_job

    :raise: RuntimeError if the validation component fails
    """
    with dsl.ParallelFor(items=validation_configs, name="validate-ingest-files-workers", parallelism=64) as item:
        components.validate_anndata_files_job(gcs_config_path=item)
