import logging
import uuid
from contextlib import contextmanager
from datetime import datetime
from typing import Any, Generator, List, Tuple

from google.api_core.exceptions import NotFound
from google.cloud import bigquery
from tenacity import before_log, retry, stop_after_attempt, wait_exponential

from cellarium.nexus.omics_datastore.bq_avro_schemas import cell_management, converter
from cellarium.nexus.omics_datastore.bq_ops import constants, exceptions
from cellarium.nexus.omics_datastore.bq_ops.create_bq_tables import create_staging_table

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def initialize_bigquery_client(project_id: str, dataset: str) -> bigquery.Client:
    """
    Initialize and return an authenticated BigQuery client.

    :param project_id: The ID of the Google Cloud project
    :param dataset: The BigQuery dataset name

    :return: An authenticated BigQuery client
    """
    client = bigquery.Client(project=project_id)
    return client


def generate_staging_suffix() -> str:
    """
    Generate a unique suffix for staging tables based on timestamp and UUID.

    :return: A unique staging suffix in format 'staging_YYYYMMDDHHMMSS_uuid'
    """
    timestamp = datetime.utcnow().strftime("%Y%m%d%H%M%S")
    unique_id = uuid.uuid4().hex[:6]
    return f"staging_{timestamp}_{unique_id}"


def define_ingestion_specs(staging_suffix: str) -> List[Tuple[str, str, str, bigquery.SourceFormat]]:
    """
    Define the specifications for each table to be ingested.

    :param staging_suffix: The unique suffix for staging tables

    :return: List of tuples containing (final_table_name, staging_table_name, source_file_pattern, file_format)
    """
    return [
        (
            constants.BQ_INGEST_TABLE_NAME,
            f"{constants.BQ_INGEST_TABLE_NAME}_{staging_suffix}",
            constants.INGEST_INGEST_FILE_NAME,
            bigquery.SourceFormat.AVRO,
        ),
        (
            constants.BQ_CELL_INFO_TABLE_NAME,
            f"{constants.BQ_CELL_INFO_TABLE_NAME}_{staging_suffix}",
            constants.INGEST_CELL_INFO_FILE_NAME,
            bigquery.SourceFormat.AVRO,
        ),
        (
            constants.BQ_FEATURE_INFO_TABLE_NAME,
            f"{constants.BQ_FEATURE_INFO_TABLE_NAME}_{staging_suffix}",
            constants.INGEST_FEATURE_INFO_FILE_NAME,
            bigquery.SourceFormat.AVRO,
        ),
        (
            constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME,
            f"{constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME}_{staging_suffix}",
            constants.INGEST_RAW_COUNTS_FILE_PATTERN,
            bigquery.SourceFormat.CSV,
        ),
    ]


def load_table_from_gcs(
    client: bigquery.Client,
    project: str,
    dataset: str,
    table_name: str,
    file_pattern: str,
    file_format: bigquery.SourceFormat,
    gcs_bucket_name: str,
    gcs_stage_dir: str,
) -> None:
    """
    Load data from GCS into a BigQuery table.

    :param client: The BigQuery client
    :param project: GCP project ID
    :param dataset: BigQuery dataset name
    :param table_name: Name of the target table
    :param file_pattern: Pattern to match source files in GCS
    :param file_format: Format of the source files (AVRO or CSV)
    :param gcs_bucket_name: Name of the GCS bucket
    :param gcs_stage_dir: Directory in GCS containing the data

    :raise Exception: If the load job fails
    """
    table_id = f"{project}.{dataset}.{table_name}"
    uri = f"gs://{gcs_bucket_name}/{gcs_stage_dir}/{file_pattern}"

    # Configure the load job
    job_config = bigquery.LoadJobConfig(source_format=file_format)

    # Start the load job
    load_job = client.load_table_from_uri(uri, table_id, job_config=job_config)
    logger.info(f"Starting job to load data from {uri} into {table_id}")

    # Wait for the job to complete
    load_job.result()
    logger.info(f"Loaded {load_job.output_rows} rows into {table_id}")


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=5, max=30),
    before=before_log(logger, logging.INFO),
)
def perform_load_table_from_gcs(
    client: bigquery.Client,
    project: str,
    dataset: str,
    table_name: str,
    file_pattern: str,
    file_format: bigquery.SourceFormat,
    gcs_bucket_name: str,
    gcs_stage_dir: str,
) -> None:
    """
    Load data from GCS into a BigQuery table with retry logic.

    :param client: The BigQuery client
    :param project: GCP project ID
    :param dataset: BigQuery dataset name
    :param table_name: Name of the target table
    :param file_pattern: Pattern to match source files in GCS
    :param file_format: Format of the source files (AVRO or CSV)
    :param gcs_bucket_name: Name of the GCS bucket
    :param gcs_stage_dir: Directory in GCS containing the data

    :raise Exception: If the load job fails after retry attempts are exhausted
    """
    # Call the original load function with all passed parameters
    load_table_from_gcs(
        client=client,
        project=project,
        dataset=dataset,
        table_name=table_name,
        file_pattern=file_pattern,
        file_format=file_format,
        gcs_bucket_name=gcs_bucket_name,
        gcs_stage_dir=gcs_stage_dir,
    )


def load_data_into_staging(
    client: bigquery.Client,
    project_id: str,
    dataset: str,
    gcs_bucket_name: str,
    gcs_stage_dir: str,
    ingestion_specs: List[Tuple[str, str, str, bigquery.SourceFormat]],
) -> bool:
    """
    Load data from GCS into staging tables with retry logic.

    :param client: The BigQuery client
    :param project_id: GCP project ID
    :param dataset: BigQuery dataset name
    :param gcs_bucket_name: GCS bucket name
    :param gcs_stage_dir: GCS staging directory
    :param ingestion_specs: List of ingestion specifications

    :return: True if all loads succeed, False otherwise
    """
    all_loads_succeeded = True

    for final_table, staging_table, file_pattern, source_format in ingestion_specs:
        try:
            logger.info(f"Loading {file_pattern} into staging table '{staging_table}'")

            perform_load_table_from_gcs(
                client=client,
                project=project_id,
                dataset=dataset,
                table_name=staging_table,
                file_pattern=file_pattern,
                file_format=source_format,
                gcs_bucket_name=gcs_bucket_name,
                gcs_stage_dir=gcs_stage_dir,
            )

            logger.info(f"Successfully loaded staging table '{staging_table}'")

        except Exception as e:
            logger.error(f"An error occurred loading into staging table {staging_table}: {e}")
            all_loads_succeeded = False

    return all_loads_succeeded


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=5, max=30),
    before=before_log(logger, logging.INFO),
)
def commit_to_production(
    client: bigquery.Client,
    dataset: str,
    ingestion_specs: List[Tuple[str, str, str, bigquery.SourceFormat]],
) -> bool:
    """
    Commit data from staging tables to production tables atomically.

    :param client: The BigQuery client
    :param dataset: BigQuery dataset name
    :param ingestion_specs: List of ingestion specifications

    :return: True if commit succeeds, False otherwise
    """
    try:
        logger.info("Beginning final transaction to move data into production tables")

        # Prepare transaction statements
        statements = ["BEGIN TRANSACTION;"]
        for final_table, staging_table, _, _ in ingestion_specs:
            statements.append(
                f"""
                INSERT INTO `{dataset}.{final_table}`
                SELECT * FROM `{dataset}.{staging_table}`;
                """
            )
        statements.append("COMMIT TRANSACTION;")

        # Execute the multi-statement transaction
        multi_statement_query = "\n".join(statements)
        logger.info("Executing multi-statement transaction for final merge...")
        job = client.query(multi_statement_query)
        job.result()  # Wait for completion

        logger.info("Transaction committed successfully. Production tables now updated.")
        return True

    except Exception as e:
        logger.error(f"Transaction failed, no data committed: {e}")
        return False


def cleanup_staging_tables(
    client: bigquery.Client,
    project_id: str,
    dataset: str,
    ingestion_specs: list[tuple[str, str, str, bigquery.SourceFormat]],
) -> None:
    """
    Delete staging tables to maintain a clean environment.

    :param client: The BigQuery client
    :param project_id: The GCP project ID
    :param dataset: The BigQuery dataset name
    :param ingestion_specs: List of ingestion specifications
    """
    logger.info("Cleaning up staging tables...")

    for _, staging_table, _, _ in ingestion_specs:
        staging_table_id = f"{project_id}.{dataset}.{staging_table}"
        try:
            client.delete_table(staging_table_id)
            logger.info(f"Dropped staging table '{staging_table_id}'")
        except NotFound:
            logger.warning(f"Staging table '{staging_table_id}' does not exist")
        except Exception as e:
            logger.warning(f"Failed to drop staging table '{staging_table_id}': {e}")


def get_schema_for_table(final_table: str) -> list[bigquery.SchemaField]:
    schema_map = {
        constants.BQ_INGEST_TABLE_NAME: cell_management.IngestInfoBQAvroSchema,
        constants.BQ_CELL_INFO_TABLE_NAME: cell_management.CellInfoBQAvroSchema,
        constants.BQ_FEATURE_INFO_TABLE_NAME: cell_management.FeatureInfoBQAvroSchema,
        constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME: cell_management.RawCountMatrixCOOBQAvroSchema,
    }

    if final_table not in schema_map:
        raise exceptions.BigQuerySchemaError(f"No schema defined for table '{final_table}'")

    return converter.pydantic_to_bigquery(pydantic_model=schema_map[final_table])


@contextmanager
def bigquery_ingest_context(
    *,
    project_id: str,
    dataset: str,
    gcs_bucket_name: str,
    gcs_stage_dir: str,
) -> Generator[None, Any, None]:
    """
    Create a context to manage the BigQuery data ingestion process.

    Handles the full lifecycle of ingesting data from GCS to BigQuery including:
    creating staging tables, loading data, committing to production tables, and cleanup.
    Data is loaded into staging tables during context entry, but commits to main
    production tables happen only when exiting the context manager successfully.

    :param project_id: GCP project ID
    :param dataset: BigQuery dataset name
    :param gcs_bucket_name: Name of the GCS bucket containing data files
    :param gcs_stage_dir: Directory within GCS bucket containing the data files

    :raise exceptions.BigQueryLoadError: If loading data to staging tables fails
    :raise exceptions.BigQueryCommitError: If committing data to production fails
    :raise Exception: If any other error occurs during the ingestion process

    :yield: None
    """
    client = initialize_bigquery_client(project_id, dataset)
    staging_suffix = generate_staging_suffix()
    ingestion_specs = define_ingestion_specs(staging_suffix)

    try:
        logger.info("Creating staging tables...")
        for final_table, _, _, _ in ingestion_specs:
            schema = get_schema_for_table(final_table)
            clustering_fields = ["cell_id"] if final_table == constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME else None

            create_staging_table(
                client=client,
                project=project_id,
                dataset=dataset,
                base_table_name=final_table,
                staging_suffix=staging_suffix,
                schema=schema,
                clustering_fields=clustering_fields,
            )

        logger.info("Loading data into staging tables...")
        all_loads_succeeded = load_data_into_staging(
            client=client,
            project_id=project_id,
            dataset=dataset,
            gcs_bucket_name=gcs_bucket_name,
            gcs_stage_dir=gcs_stage_dir,
            ingestion_specs=ingestion_specs,
        )

        if not all_loads_succeeded:
            raise exceptions.BigQueryLoadError("Not all staging loads succeeded. Skipping commit to production.")

        yield

        logger.info("Committing data to production tables...")
        commit_success = commit_to_production(
            client=client,
            dataset=dataset,
            ingestion_specs=ingestion_specs,
        )

        if not commit_success:
            raise exceptions.BigQueryCommitError("Commit to production failed. No data was appended.")

        logger.info("BigQuery ingest committed to production successfully.")

    except Exception as e:
        logger.error(f"BigQuery ingest failed: {e}")
        raise

    finally:
        try:
            cleanup_staging_tables(
                client=client,
                project_id=project_id,
                dataset=dataset,
                ingestion_specs=ingestion_specs,
            )
        except Exception as cleanup_error:
            logger.error(f"Failed to cleanup staging tables: {cleanup_error}")
