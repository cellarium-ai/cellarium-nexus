"""
Prepare BigQuery tables for efficient data extraction by preprocessing and staging data.
"""

import csv
import logging
import tempfile
from datetime import datetime, timedelta, timezone

from pathlib import Path
from typing import Any, Sequence

from google.cloud import bigquery

from cellarium.nexus.omics_datastore import bq_sql
from cellarium.nexus.omics_datastore.bq_ops import constants
from cellarium.nexus.omics_datastore.bq_ops.extract.metadata_extractor import MetadataExtractor
from cellarium.nexus.shared import schemas

logger = logging.getLogger(__name__)

# Template paths
TEMPLATE_DIR = Path(__file__).parent.parent.parent / "sql_templates" / "prepare_extract"
CELL_INFO_RAND_TEMPLATE = TEMPLATE_DIR / "prepare_cell_info_randomized.sql.mako"
CELL_INFO_TEMPLATE = TEMPLATE_DIR / "prepare_cell_info.sql.mako"
DROP_CELL_INFO_RAND_TEMPLATE = TEMPLATE_DIR / "drop_prepare_cell_info_randomized.sql.mako"


class ExtractTablePreparer:
    """
    Prepares BigQuery tables for efficient data extraction.
    """

    def __init__(
        self,
        client: bigquery.Client,
        project: str,
        dataset: str,
        extract_table_prefix: str,
    ):
        """
        Initialize extract table preparer.

        :param client: BigQuery client instance
        :param project: GCP project ID
        :param dataset: BigQuery dataset ID
        :param extract_table_prefix: Prefix for extract table names
        """
        self.client = client
        self.project = project
        self.dataset = dataset
        self.prefix = extract_table_prefix

    def execute_query(self, sql: str) -> Any:
        """
        Execute a BigQuery SQL query.

        :param sql: SQL query string to execute

        :raise google.api_core.exceptions.GoogleAPIError: If query execution fails

        :return: Query results
        """
        query = self.client.query(sql)
        return query.result()

    def set_table_expiration(
        self,
        table_id: str,
        expiration_delta: timedelta,
    ) -> None:
        """
        Set expiration time for a BigQuery table.

        :param table_id: Fully qualified table ID (project.dataset.table)
        :param expiration_delta: Time delta after which the table will expire

        :raise google.api_core.exceptions.GoogleAPIError: If table update fails
        """
        expiration = datetime.now(timezone.utc) + expiration_delta
        table = self.client.get_table(table_id)
        table.expires = expiration
        self.client.update_table(table=table, fields=["expires"])
        logger.info(f"Set expiration for {table_id} to {expiration.isoformat()} UTC")

    def prepare_feature_table(self, features: Sequence[schemas.FeatureSchema]) -> None:
        """
        Create feature mapping table from provided schema using CSV load job.

        :param features: Sequence of feature schema objects

        :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
        """
        schema = [
            bigquery.SchemaField("id", "INTEGER"),
            bigquery.SchemaField("symbol", "STRING"),
            bigquery.SchemaField("ensemble_id", "STRING"),
        ]

        # Create temporary CSV file
        with tempfile.NamedTemporaryFile(mode="w", newline="", suffix=".csv") as temp_file:
            writer = csv.writer(temp_file)
            # Write features to CSV
            for feature in features:
                writer.writerow([feature.id, feature.symbol, feature.ensemble_id])
            temp_file.flush()

            # Configure load job
            job_config = bigquery.LoadJobConfig(
                schema=schema,
                skip_leading_rows=0,
                source_format=bigquery.SourceFormat.CSV,
                write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
            )

            table_id = f"{self.project}.{self.dataset}.{self.prefix}{constants.BQ_EXTRACT_FEATURE_INFO_TABLE_NAME}"

            # Load data from CSV
            logger.info("Creating feature info table...")
            with open(temp_file.name, "rb") as source_file:
                load_job = self.client.load_table_from_file(source_file, table_id, job_config=job_config)
                load_job.result()  # Wait for job to complete

            logger.info(f"Loaded {load_job.output_rows} rows into {table_id}")

            # Set the expiration time
            self.set_table_expiration(table_id=table_id, expiration_delta=timedelta(hours=3))

    def prepare_cell_info(
        self,
        extract_bin_size: int | None = None,
        random_seed_offset: int = 0,
        partition_bin_count: int = 40000,
        partition_size: int = 10,
        extract_bin_keys: list[str] | None = None,
        filters: dict[str, Any] | None = None,
        obs_columns: list[str] | None = None,
        metadata_extra_columns: list[str] | None = None,
    ) -> None:
        """
        Create cell info table with binning and randomization.

        :param extract_bin_size: Size of cell bins
        :param random_seed_offset: Offset for randomization
        :param partition_bin_count: Number of partitions
        :param partition_size: Size of each partition
        :param extract_bin_keys: List of keys to bin by. If not provided, bins will be assigned randomly.
        :param filters: Query filters to apply
        :param obs_columns: Observation columns to include
        :param metadata_extra_columns: Additional metadata columns to include to extract files from `metadata_extra`
            JSON blob. If not provided, none will be included.

        :raise ValueError: If binning parameters are invalid
        :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
        """
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
            select=obs_columns,
            filters=filters,
            random_seed_offset=random_seed_offset,
            metadata_columns=metadata_extra_columns,
            extract_bin_keys=extract_bin_keys,
        )

        # Create randomized table
        logger.info("Creating randomized cell info table...")
        sql = bq_sql.render(str(CELL_INFO_RAND_TEMPLATE), template_data)
        self.execute_query(sql)

        # Create binned table
        logger.info("Creating binned cell info table...")
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
            select=obs_columns,
            filters=filters,
            random_seed_offset=random_seed_offset,
            partition_bin_count=partition_bin_count,
            partition_size=partition_size,
            extract_bin_size=extract_bin_size,
            extract_bin_keys=extract_bin_keys,
            metadata_columns=metadata_extra_columns,
        )
        sql = bq_sql.render(str(CELL_INFO_TEMPLATE), template_data)
        self.execute_query(sql)

        # Set expiration time for the cell info table (14 days)
        cell_info_table_id = f"{self.project}.{self.dataset}.{self.prefix}{constants.BQ_EXTRACT_CELL_INFO_TABLE_NAME}"
        self.set_table_expiration(table_id=cell_info_table_id, expiration_delta=timedelta(days=14))

        # Drop intermediate table
        logger.info("Dropping intermediate randomized table...")
        sql = bq_sql.render(str(DROP_CELL_INFO_RAND_TEMPLATE), template_data)
        self.execute_query(sql)

    def prepare_count_matrix(self, partition_bin_count: int = 40000, partition_size: int = 10) -> None:
        """
        Create count matrix table with feature data arrays.

        :param partition_bin_count: Number of partitions
        :param partition_size: Size of each partition

        :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
        """
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
            partition_bin_count=partition_bin_count,
            partition_size=partition_size,
        )

        logger.info("Creating count matrix table...")
        sql = bq_sql.render(str(TEMPLATE_DIR / "prepare_count_matrix.sql.mako"), template_data)
        self.execute_query(sql)

        # Set expiration time for the count matrix table (4 hours)
        count_matrix_table_id = (
            f"{self.project}.{self.dataset}.{self.prefix}{constants.BQ_EXTRACT_MATRIX_COO_TABLE_NAME}"
        )
        self.set_table_expiration(table_id=count_matrix_table_id, expiration_delta=timedelta(hours=4))


def prepare_extract_tables(
    client: bigquery.Client,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    features: Sequence[schemas.FeatureSchema],
    categorical_column_count_limit: int = constants.CATEGORICAL_COLUMN_COUNT_LIMIT_DEFAULT,
    extract_bin_size: int | None = None,
    random_seed_offset: int = 0,
    partition_bin_count: int = 40000,
    partition_size: int = 10,
    extract_bin_keys: list[str] | None = None,
    filters: dict[str, Any] | None = None,
    obs_columns: list[str] | None = None,
    metadata_extra_columns: list[str] | None = None,
) -> schemas.ExtractMetadata:
    """
    Prepare all necessary tables for data extraction.

    :param client: BigQuery client instance
    :param project: GCP project ID
    :param dataset: BigQuery dataset ID
    :param extract_table_prefix: Prefix for extract table names
    :param features: Sequence of feature schema objects
    :param categorical_column_count_limit: Maximum number of categories per categorical column to be considered as
        categorical. If the number of categories exceeds this limit, the column will not be unified across all extract
        files.
    :param extract_bin_size: Size of cell bins
    :param random_seed_offset: Offset for randomization
    :param partition_bin_count: Number of partitions
    :param partition_size: Size of each partition
    :param extract_bin_keys: List of keys to bin by. If not provided, bins will be assigned randomly.
    :param filters: Query filters to apply
    :param obs_columns: Observation columns to include
    :param metadata_extra_columns: Additional metadata columns to include to extract files from `metadata_extra`.
        If not provided, none will be included.

    :raise ValueError: If binning parameters are invalid
    :raise google.api_core.exceptions.GoogleAPIError: If table creation fails

    :return: The complete ExtractMetadata object
    """
    preparer = ExtractTablePreparer(
        client=client,
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
    )

    logger.info("Starting extract table preparation...")

    # Step 1: Prepare feature mapping table
    logger.info("Preparing feature mapping table...")
    preparer.prepare_feature_table(features)

    # Step 2: Prepare cell info with binning
    logger.info("Preparing cell info table with binning...")
    preparer.prepare_cell_info(
        extract_bin_size=extract_bin_size,
        random_seed_offset=random_seed_offset,
        partition_bin_count=partition_bin_count,
        partition_size=partition_size,
        extract_bin_keys=extract_bin_keys,
        filters=filters,
        obs_columns=obs_columns,
        metadata_extra_columns=metadata_extra_columns,
    )

    # Step 3: Prepare count matrix
    logger.info("Preparing count matrix table...")
    preparer.prepare_count_matrix(
        partition_bin_count=partition_bin_count,
        partition_size=partition_size,
    )
    logger.info("Extract table preparation completed successfully")

    # Step 4: Extract metadata
    logger.info("Extracting metadata...")
    metadata_extractor = MetadataExtractor(
        client=client,
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
        filters=filters,
        extract_bin_size=extract_bin_size,
        categorical_column_count_limit=categorical_column_count_limit,
    )

    # Compose the extract metadata
    return metadata_extractor.compose_extract_metadata()
