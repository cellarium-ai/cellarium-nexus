import logging

from google.api_core.exceptions import Conflict
from google.cloud import bigquery
from nexus.omics_datastore.bq_avro_schemas import cell_management, converter
from nexus.omics_datastore.bq_ops import constants

# Configure logging
logger = logging.getLogger(__name__)


def create_table(
    client: bigquery.Client,
    project: str,
    dataset: str,
    table_name: str,
    schema: list[bigquery.SchemaField],
    clustering_fields: list[str] | None = None,
) -> None:
    """
    Create the specified table in the specified project / dataset using the specified schema and clustering fields.

    :param client: BigQuery client instance
    :param project: GCP project ID
    :param dataset: BigQuery dataset name
    :param table_name: Name of the table to create
    :param schema: List of SchemaField objects defining the table schema
    :param clustering_fields: Optional list of field names to use for clustering

    :raise Conflict: if table already exists
    :raise Exception: for other BigQuery errors
    """
    table_id = f"{project}.{dataset}.{table_name}"

    # Create table object with specified schema
    table = bigquery.Table(table_ref=table_id, schema=schema)

    # Add clustering if specified
    if clustering_fields is not None:
        table.clustering_fields = clustering_fields

    try:
        _ = client.create_table(table)
        logger.info(f"Created '{table_id}'.")
    except Conflict:
        logger.info(f"Table '{table_id}' exists, continuing.")


def create_dataset(client: bigquery.Client, project: str, dataset: str, location: str) -> None:
    """
    Create the specified dataset in the specified project and location.

    :param client: BigQuery client instance
    :param project: GCP project ID
    :param dataset: Name of the dataset to create
    :param location: Geographic location for the dataset (e.g., "US")

    :raise Conflict: if dataset already exists
    :raise Exception: for other BigQuery errors
    """
    dataset_id = f"{project}.{dataset}"

    # Construct a full Dataset object to send to the API
    dataset = bigquery.Dataset(dataset_ref=dataset_id)
    dataset.location = location

    try:
        # Make an API request with an explicit timeout
        _ = client.create_dataset(dataset, timeout=30)
        logger.info(f"Created dataset {dataset_id}.")
    except Conflict:
        logger.info(f"Dataset {dataset_id} exists, continuing.")


def create_bigquery_objects(client: bigquery.Client, project: str, dataset: str, location: str) -> str:
    """
    Create the core BigQuery dataset and tables for cell management.

    :param client: BigQuery client instance
    :param project: GCP project ID
    :param dataset: Name of the dataset to create
    :param location: Geographic location for the dataset (e.g., "US")

    :raise Exception: if any BigQuery operation fails

    :return: URL to the created dataset in BigQuery console
    """
    # Create the dataset first
    create_dataset(client=client, project=project, dataset=dataset, location=location)

    # Convert Pydantic models to BigQuery schemas
    bq_ingest_info_table_schema = converter.pydantic_to_bigquery(pydantic_model=cell_management.IngestInfoBQAvroSchema)
    bq_cell_info_table_schema = converter.pydantic_to_bigquery(pydantic_model=cell_management.CellInfoBQAvroSchema)
    bq_feature_info_table_schema = converter.pydantic_to_bigquery(
        pydantic_model=cell_management.FeatureInfoBQAvroSchema
    )
    bq_raw_matrix_table_schema = converter.pydantic_to_bigquery(
        pydantic_model=cell_management.RawCountMatrixCOOBQAvroSchema
    )

    # Create all required tables
    create_table(
        client=client,
        project=project,
        dataset=dataset,
        table_name=constants.BQ_INGEST_TABLE_NAME,
        schema=bq_ingest_info_table_schema,
    )
    create_table(
        client=client,
        project=project,
        dataset=dataset,
        table_name=constants.BQ_CELL_INFO_TABLE_NAME,
        schema=bq_cell_info_table_schema,
    )
    create_table(
        client=client,
        project=project,
        dataset=dataset,
        table_name=constants.BQ_FEATURE_INFO_TABLE_NAME,
        schema=bq_feature_info_table_schema,
    )
    create_table(
        client=client,
        project=project,
        dataset=dataset,
        table_name=constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME,
        schema=bq_raw_matrix_table_schema,
        clustering_fields=["cell_id"],  # Optimize queries that filter by cell_id
    )

    # Return URL to the dataset in BigQuery console
    return f"https://console.cloud.google.com/bigquery?project={project}&p={project}&d={dataset}&page=dataset"


def create_staging_table(
    client: bigquery.Client,
    project: str,
    dataset: str,
    base_table_name: str,
    staging_suffix: str,
    schema: list[bigquery.SchemaField],
    clustering_fields: list[str] | None = None,
) -> None:
    """
    Create a staging table within the specified dataset using a unique suffix to avoid conflicts.

    :param client: BigQuery client instance
    :param project: GCP project ID
    :param dataset: BigQuery dataset name
    :param base_table_name: The base name of the table to create
    :param staging_suffix: A unique suffix to append to the table name
    :param schema: List of SchemaField objects defining the table schema
    :param clustering_fields: Optional list of field names to use for clustering

    :raise Conflict: if table already exists
    :raise Exception: for other BigQuery errors
    """
    staging_table_name = f"{base_table_name}_{staging_suffix}"
    table_id = f"{project}.{dataset}.{staging_table_name}"

    # Create table object with specified schema
    table = bigquery.Table(table_id, schema=schema)
    if clustering_fields:
        table.clustering_fields = clustering_fields

    try:
        client.create_table(table)
        logger.info(f"Created staging table '{table_id}'.")
    except Conflict:
        logger.info(f"Staging table '{table_id}' already exists, continuing.")
    except Exception as e:
        logger.error(f"Failed to create staging table '{table_id}': {e}")
        raise
