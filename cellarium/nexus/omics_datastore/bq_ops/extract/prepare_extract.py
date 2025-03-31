"""
Prepare BigQuery tables for efficient data extraction by preprocessing and staging data.
"""

import csv
import logging
import tempfile
from pathlib import Path
from typing import Any, Sequence

from google.cloud import bigquery
from pydantic import BaseModel
from nexus.omics_datastore import bq_sql
from nexus.omics_datastore.bq_ops.extract.metadata_extractor import MetadataExtractor

logger = logging.getLogger(__name__)

# Template paths
TEMPLATE_DIR = Path(__file__).parent.parent.parent / "sql_templates" / "prepare_extract"
CELL_INFO_RAND_TEMPLATE = TEMPLATE_DIR / "prepare_cell_info_randomized.sql.mako"
CELL_INFO_TEMPLATE = TEMPLATE_DIR / "prepare_cell_info.sql.mako"
DROP_CELL_INFO_RAND_TEMPLATE = TEMPLATE_DIR / "drop_prepare_cell_info_randomized.sql.mako"


class FeatureSchema(BaseModel):
    """
    Schema for feature data.

    :param id: Unique identifier for the feature
    :param symbol: Gene symbol
    :param ensemble_id: Ensemble identifier
    """

    id: int
    symbol: str
    ensemble_id: str


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

    def prepare_feature_table(self, features: Sequence[FeatureSchema]) -> None:
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

            table_id = f"{self.project}.{self.dataset}.{self.prefix}__extract_feature_info"

            # Load data from CSV
            logger.info("Creating feature info table...")
            with open(temp_file.name, "rb") as source_file:
                load_job = self.client.load_table_from_file(source_file, table_id, job_config=job_config)
                load_job.result()  # Wait for job to complete

            logger.info(f"Loaded {load_job.output_rows} rows into {table_id}")

    def prepare_cell_info(
        self,
        extract_bin_size: int | None = None,
        assign_bin_by_category: bool = False,
        extract_bin_category_column_name: str | None = None,
        random_seed_offset: int = 0,
        partition_bin_count: int = 40000,
        partition_size: int = 10,
        filters: dict[str, Any] | None = None,
        obs_columns: list[str] | None = None,
    ) -> None:
        """
        Create cell info table with binning and randomization.

        :param extract_bin_size: Size of cell bins
        :param assign_bin_by_category: Whether to bin by category
        :param extract_bin_category_column_name: Column name for category binning
        :param random_seed_offset: Offset for randomization
        :param partition_bin_count: Number of partitions
        :param partition_size: Size of each partition
        :param filters: Query filters to apply
        :param obs_columns: Observation columns to include

        :raise ValueError: If binning parameters are invalid
        :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
        """
        if not assign_bin_by_category and extract_bin_size is None:
            raise ValueError("extract_bin_size required when not binning by category")
        if assign_bin_by_category and extract_bin_category_column_name is None:
            raise ValueError("category column name required when binning by category")

        # Create randomized intermediate table
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
            select=obs_columns,
            filters=filters,
            random_seed_offset=random_seed_offset,
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
            assign_bin_by_category=assign_bin_by_category,
            extract_bin_category_column_name=extract_bin_category_column_name,
        )
        sql = bq_sql.render(str(CELL_INFO_TEMPLATE), template_data)
        self.execute_query(sql)

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


def prepare_extract_tables(
    client: bigquery.Client,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    features: Sequence[FeatureSchema],
    extract_bin_size: int | None = None,
    assign_bin_by_category: bool = False,
    extract_bin_category_column_name: str | None = None,
    random_seed_offset: int = 0,
    partition_bin_count: int = 40000,
    partition_size: int = 10,
    filters: dict[str, Any] | None = None,
    obs_columns: list[str] | None = None,
    bucket_name: str | None = None,
    extract_bucket_path: str | None = None,
) -> None:
    """
    Prepare all necessary tables for data extraction.

    :param client: BigQuery client instance
    :param project: GCP project ID
    :param dataset: BigQuery dataset ID
    :param extract_table_prefix: Prefix for extract table names
    :param features: Sequence of feature schema objects
    :param extract_bin_size: Size of cell bins
    :param assign_bin_by_category: Whether to bin by category
    :param extract_bin_category_column_name: Column name for category binning
    :param random_seed_offset: Offset for randomization
    :param partition_bin_count: Number of partitions
    :param partition_size: Size of each partition
    :param filters: Query filters to apply
    :param obs_columns: Observation columns to include
    :param bucket_name: GCS bucket name for metadata storage
    :param extract_bucket_path: Path within bucket for metadata storage

    :raise ValueError: If binning parameters are invalid
    :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
    :raise IOError: If metadata file operations fail
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
        assign_bin_by_category=assign_bin_by_category,
        extract_bin_category_column_name=extract_bin_category_column_name,
        random_seed_offset=random_seed_offset,
        partition_bin_count=partition_bin_count,
        partition_size=partition_size,
        filters=filters,
        obs_columns=obs_columns,
    )

    # Step 3: Prepare count matrix
    logger.info("Preparing count matrix table...")
    preparer.prepare_count_matrix(
        partition_bin_count=partition_bin_count,
        partition_size=partition_size,
    )
    logger.info("Extract table preparation completed successfully")

    # Step 4: Extract and save metadata if storage location provided
    if bucket_name and extract_bucket_path:
        logger.info("Extracting and saving metadata...")
        metadata_extractor = MetadataExtractor(
            client=client,
            project=project,
            dataset=dataset,
            extract_table_prefix=extract_table_prefix,
            filters=filters,
        )
        metadata_extractor.save_metadata(
            bucket_name=bucket_name,
            extract_bucket_path=extract_bucket_path,
        )
        logger.info("Metadata extraction completed successfully")
