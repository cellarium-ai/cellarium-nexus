"""
Provide a centralized controller for BigQuery operations in the Nexus datastore.
"""

import logging
from pathlib import Path
from typing import Any, ContextManager, Sequence

from google.cloud import bigquery

from cellarium.nexus.omics_datastore.bq_ops import bq_sql, constants
from cellarium.nexus.omics_datastore.bq_ops.create_bq_tables import create_bigquery_objects
from cellarium.nexus.omics_datastore.bq_ops.extract.extract import extract_bins
from cellarium.nexus.omics_datastore.bq_ops.extract.metadata_extractor import MetadataExtractor
from cellarium.nexus.omics_datastore.bq_ops.extract.prepare_extract import prepare_extract_tables
from cellarium.nexus.omics_datastore.bq_ops.ingest.create_ingest_files import create_ingest_files
from cellarium.nexus.omics_datastore.bq_ops.ingest.ingest_data_to_bigquery import bigquery_ingest_context
from cellarium.nexus.shared import schemas

logger = logging.getLogger(__name__)


class BigQueryDataOperator:
    """
    Control and manage BigQuery operations for the Nexus datastore.

    This class provides low-level operations for interacting with BigQuery:
    - Dataset and table management
    - Data ingestion
    - Extract table preparation
    - Metadata extraction
    - Data extraction
    """

    def __init__(
        self,
        *,
        client: bigquery.Client,
        project: str,
        dataset: str,
    ) -> None:
        """
        Initialize the BigQuery datastore operator.

        :param client: Initialized BigQuery client
        :param project: GCP project ID
        :param dataset: BigQuery dataset name

        :raise ValueError: If any of the required parameters are invalid
        """
        self.client = client
        self.project = project
        self.dataset = dataset

    def initialize_bigquery_resources(self, *, location: str = "US", labels: dict[str, str] | None = None) -> str:
        """
        Create BigQuery dataset and all required tables.

        :param location: GCP region where to create the resources
        :param labels: Optional dictionary of labels to apply to the dataset and tables

        :raise google.api_core.exceptions.GoogleAPIError: If resource creation fails
        """
        return create_bigquery_objects(
            client=self.client, project=self.project, dataset=self.dataset, location=location, labels=labels
        )

    def create_ingest_files(
        self,
        *,
        adata_file_path: Path,
        tag: str | None,
        cell_info_start_index: int,
        cell_info_end_index: int,
        feature_info_start_index: int,
        feature_info_end_index: int,
        ingest_id: int,
        output_dir: Path,
        column_mapping: dict[str, Any] | None = None,
        uns_keys_to_keep: list[str] | None = None,
    ) -> dict[str, Any]:
        """
        Create ingest files for BigQuery ingestion.

        :param adata_file_path: Path to the AnnData file
        :param tag: Optional tag for the ingest
        :param cell_info_start_index: Start index for cell info
        :param cell_info_end_index: End index for cell info
        :param feature_info_start_index: Start index for feature info
        :param feature_info_end_index: End index for feature info
        :param ingest_id: ID of the ingest
        :param output_dir: Directory to write output files
        :param column_mapping: Optional dictionary containing obs and var column mappings
        :param uns_keys_to_keep: Optional list of keys to keep in the `uns` JSON blob. If not provided,
            all keys will be kept

        :raise Exception: If file creation fails

        :return: Dictionary containing ingest job results
        """
        return create_ingest_files(
            adata_file_path=adata_file_path,
            tag=tag,
            cell_info_start_index=cell_info_start_index,
            cell_info_end_index=cell_info_end_index,
            feature_info_start_index=feature_info_start_index,
            feature_info_end_index=feature_info_end_index,
            ingest_id=ingest_id,
            output_dir=output_dir,
            column_mapping=column_mapping,
            uns_keys_to_keep=uns_keys_to_keep,
        )

    def ingest_data(
        self,
        *,
        gcs_bucket_name: str,
        gcs_stage_dir: str,
    ) -> ContextManager:
        """
        Ingest data from GCS into BigQuery tables.

        :param gcs_bucket_name: GCS bucket name containing the data
        :param gcs_stage_dir: Directory in the bucket containing staged files

        :raise google.api_core.exceptions.GoogleAPIError: If ingestion fails

        :return: Context manager for the BigQuery ingestion job.
        """
        return bigquery_ingest_context(
            project_id=self.project, dataset=self.dataset, gcs_bucket_name=gcs_bucket_name, gcs_stage_dir=gcs_stage_dir
        )

    def extract_metadata(
        self,
        *,
        extract_table_prefix: str,
        extract_bin_size: int,
        filters: dict[str, Any] | None = None,
    ) -> MetadataExtractor:
        """
        Create a metadata extractor for the given extract tables.

        :param extract_table_prefix: Prefix for extract table names
        :param extract_bin_size: Size for extract bins
        :param filters: Optional filters to apply during extraction

        :return: Initialized metadata extractor
        """
        return MetadataExtractor(
            client=self.client,
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=extract_table_prefix,
            filters=filters,
            extract_bin_size=extract_bin_size,
        )

    def prepare_extract_tables(
        self,
        *,
        extract_table_prefix: str,
        features: Sequence[schemas.FeatureSchema],
        categorical_column_count_limit: int,
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
        Prepare BigQuery tables for efficient data extraction.

        :param extract_table_prefix: Prefix for extract table names
        :param features: List of feature schemas to prepare
        :param categorical_column_count_limit: Maximum number of categories per categorical column to be considered as
            categorical. If the number of categories exceeds this limit, the column will not be unified across all
            extract files.
        :param extract_bin_size: Optional size for extract bins
        :param random_seed_offset: Offset for random seed
        :param partition_bin_count: Number of bins per partition
        :param partition_size: Size of each partition
        :param extract_bin_keys: Optional list of keys to bin by. If not provided, bins will be assigned randomly.
        :param filters: Optional filters to apply during preparation
        :param obs_columns: Optional list of observation columns to include
        :param metadata_extra_columns: Optional list of metadata extra columns to include to extract files from
            JSON blob. If not provided, none will be included.

        :raise google.api_core.exceptions.GoogleAPIError: If table preparation fails

        :return: The complete ExtractMetadata object
        """

        return prepare_extract_tables(
            client=self.client,
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=extract_table_prefix,
            features=features,
            categorical_column_count_limit=categorical_column_count_limit,
            extract_bin_size=extract_bin_size,
            random_seed_offset=random_seed_offset,
            partition_bin_count=partition_bin_count,
            partition_size=partition_size,
            extract_bin_keys=extract_bin_keys,
            filters=filters,
            obs_columns=obs_columns,
            metadata_extra_columns=metadata_extra_columns,
        )

    def extract_data(
        self,
        *,
        extract_table_prefix: str,
        bins: list[int],
        output_dir: Path,
        extract_metadata: schemas.ExtractMetadata,
        obs_columns: list[str] | None = None,
        max_workers: int | None = None,
    ) -> None:
        """
        Extract data from prepared extract tables into AnnData files.

        :param extract_table_prefix: Prefix for extract table names
        :param bins: List of bin numbers to extract
        :param output_dir: Local directory to save AnnData files
        :param extract_metadata: ExtractMetadata instance with metadata for the extract
        :param obs_columns: Optional list of observation columns to include
        :param max_workers: Maximum number of parallel workers

        :raise google.api_core.exceptions.GoogleAPIError: If extraction fails
        :raise IOError: If file operations fail
        """
        extract_bins(
            client=self.client,
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=extract_table_prefix,
            bins=bins,
            output_dir=output_dir,
            obs_columns=obs_columns,
            max_workers=max_workers,
            extract_metadata=extract_metadata,
        )

    def count_cells(self, *, filter_statements: dict[str, Any] | None = None, dataset: str | None = None) -> int:
        """
        Count the number of cells in the cell_info table matching the given filters.

        :param filter_statements: Optional dictionary of filters to apply (column_name: value)
        :param dataset: Optional BigQuery dataset name to query (defaults to controller's dataset)

        :raise google.api_core.exceptions.GoogleAPIError: If query execution fails

        :return: The total count of matching cells
        """
        target_dataset = dataset if dataset else self.dataset
        template_path = Path(__file__).parent / "sql_templates" / "general" / "count_cells.sql.mako"

        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=target_dataset,
            table_name=constants.BQ_CELL_INFO_TABLE_NAME,
            filter_statements=filter_statements,
        )

        sql = bq_sql.render(str(template_path), template_data)
        logger.info(f"Executing cell count query on {target_dataset}...")
        query_job = self.client.query(sql)
        results = query_job.result()  # Waits for the query to finish

        # Extract the single count value
        count = 0
        for row in results:
            count = row.total_cells
            break  # Should only be one row

        logger.info(f"Found {count} matching cells.")
        return count

    def get_distinct_count(self, *, table_name: str, column_name: str, dataset: str | None = None) -> int:
        """
        Count distinct values for a column in a base table.

        :param table_name: Base table name to query
        :param column_name: Column to count distinct values for
        :param dataset: Optional override dataset; defaults to operator's dataset

        :raise google.api_core.exceptions.GoogleAPIError: If query execution fails

        :return: Distinct values count for the column
        """
        target_dataset = dataset if dataset else self.dataset
        template_path = Path(__file__).parent / "sql_templates" / "general" / "column_distinct_count.sql.mako"
        template_data = bq_sql.TemplateData(
            project=self.project, dataset=target_dataset, table_name=table_name, column_name=column_name
        )
        sql = bq_sql.render(str(template_path), template_data)
        query_job = self.client.query(sql)
        results = list(query_job.result())
        return int(results[0][0]) if results else 0

    def get_distinct_obs_values(
        self,
        *,
        column_name: str,
    ) -> list[str]:
        """
        Fetch distinct values for an obs/cell_info column.

        :param column_name: Column to retrieve values for

        :raise google.api_core.exceptions.GoogleAPIError: If query execution fails

        :return: List of distinct values
        """
        template_path = Path(__file__).parent / "sql_templates" / "general" / "column_distinct_values.sql.mako"
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            table_name=constants.BQ_CELL_INFO_TABLE_NAME,
            column_name=column_name,
        )
        sql = bq_sql.render(str(template_path), template_data)
        query_job = self.client.query(sql)
        results = list(query_job.result())
        return [str(row[0]) for row in results if row[0] is not None]

    def get_categorical_obs_columns(
        self,
        *,
        threshold: int,
        exclude: list[str] | None = None,
    ) -> set[str]:
        """
        Determine categorical string columns in obs/cell_info by distinct count threshold.

        :param threshold: Maximum distinct values to consider a column categorical
        :param exclude: Optional list of column names to skip

        :raise google.api_core.exceptions.GoogleAPIError: If schema or query fails

        :return: Set of column names that are categorical
        """
        exclude_set = set(exclude or [])
        table_ref = f"{self.project}.{self.dataset}.{constants.BQ_CELL_INFO_TABLE_NAME}"
        table = self.client.get_table(table_ref)

        # Collect string columns to evaluate
        string_columns: list[str] = [
            f.name for f in table.schema if f.field_type == "STRING" and f.name not in exclude_set
        ]

        if not string_columns:
            return set()

        # Use a single query to get distinct counts for all columns
        counts = self._get_multi_distinct_counts(
            table_name=constants.BQ_CELL_INFO_TABLE_NAME,
            column_names=string_columns,
        )

        categorical: set[str] = set()
        for col, cnt in counts.items():
            if int(cnt) <= threshold:
                categorical.add(col)
        return categorical

    def _get_multi_distinct_counts(
        self,
        *,
        table_name: str,
        column_names: list[str],
    ) -> dict[str, int]:
        """
        Count distinct values for multiple columns in a single query.

        :param table_name: Base table name to query
        :param column_names: List of column names to include in the count

        :raise google.api_core.exceptions.GoogleAPIError: If query execution fails

        :return: Mapping of column name -> distinct count
        """
        if not column_names:
            return {}

        template_path = Path(__file__).parent / "sql_templates" / "general" / "multi_column_distinct_counts.sql.mako"
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            table_name=table_name,
            column_names=column_names,
        )
        sql = bq_sql.render(str(template_path), template_data)
        query_job = self.client.query(sql)
        results = list(query_job.result())
        if not results:
            return {}

        row = results[0]
        out: dict[str, int] = {}
        # Aliases are produced as distinct_<column>
        for col in column_names:
            alias = f"distinct_{col}"
            out[col] = int(row[alias]) if alias in row else int(row.get(alias, 0))
        return out
