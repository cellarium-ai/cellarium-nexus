"""
Provide a centralized controller for BigQuery operations in the Nexus datastore.
"""

import logging
from pathlib import Path
from typing import Any, List, Optional, Sequence

from google.cloud import bigquery
from nexus.omics_datastore.bq_ops.create_bq_tables import create_bigquery_objects
from nexus.omics_datastore.bq_ops.extract.extract import extract_bins
from nexus.omics_datastore.bq_ops.extract.metadata_extractor import MetadataExtractor
from nexus.omics_datastore.bq_ops.extract.prepare_extract import FeatureSchema, prepare_extract_tables
from nexus.omics_datastore.bq_ops.ingest.create_ingest_files import create_ingest_files
from nexus.omics_datastore.bq_ops.ingest.ingest_data_to_bigquery import ingest_data_to_bigquery

logger = logging.getLogger(__name__)


class BQDatastoreController:
    """
    Control and manage BigQuery operations for the Nexus datastore.

    This controller provides low-level operations for interacting with BigQuery:
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
        Initialize the BigQuery datastore controller.

        :param client: Initialized BigQuery client
        :param project: GCP project ID
        :param dataset: BigQuery dataset name

        :raise ValueError: If any of the required parameters are invalid
        """
        self.client = client
        self.project = project
        self.dataset = dataset

    def create_dataset(self, *, location: str = "US") -> str:
        """
        Create a new BigQuery dataset.

        :param location: GCP region where to create the dataset

        :raise google.api_core.exceptions.GoogleAPIError: If dataset creation fails

        :return: Full path to the created dataset
        """
        dataset_ref = f"{self.project}.{self.dataset}"
        dataset = bigquery.Dataset(dataset_ref)
        dataset.location = location

        self.client.create_dataset(dataset, exists_ok=True)
        return dataset_ref

    def initialize_bigquery_resources(self, *, location: str = "US") -> str:
        """
        Create BigQuery dataset and all required tables.

        :param location: GCP region where to create the resources

        :raise google.api_core.exceptions.GoogleAPIError: If resource creation fails
        """
        return create_bigquery_objects(
            client=self.client,
            project=self.project,
            dataset=self.dataset,
            location=location,
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
        )

    def ingest_data(
        self,
        *,
        gcs_bucket_name: str,
        ingest_id: str,
        gcs_stage_dir: str,
    ) -> None:
        """
        Ingest data from GCS into BigQuery tables.

        :param gcs_bucket_name: GCS bucket name containing the data
        :param ingest_id: Unique identifier for this ingestion
        :param gcs_stage_dir: Directory in the bucket containing staged files

        :raise google.api_core.exceptions.GoogleAPIError: If ingestion fails
        """
        ingest_data_to_bigquery(
            project_id=self.project,
            dataset=self.dataset,
            gcs_bucket_name=gcs_bucket_name,
            ingest_id=ingest_id,
            gcs_stage_dir=gcs_stage_dir,
        )

    def extract_metadata(
        self,
        *,
        extract_table_prefix: str,
        filters: dict[str, Any] | None = None,
    ) -> MetadataExtractor:
        """
        Create a metadata extractor for the given extract tables.

        :param extract_table_prefix: Prefix for extract table names
        :param filters: Optional filters to apply during extraction

        :return: Initialized metadata extractor
        """
        return MetadataExtractor(
            client=self.client,
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=extract_table_prefix,
            filters=filters,
        )

    def prepare_extract_tables(
        self,
        *,
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
        Prepare BigQuery tables for efficient data extraction.

        :param extract_table_prefix: Prefix for extract table names
        :param features: List of feature schemas to prepare
        :param extract_bin_size: Optional size for extract bins
        :param assign_bin_by_category: Whether to assign bins by category
        :param extract_bin_category_column_name: Column name for category-based binning
        :param random_seed_offset: Offset for random seed
        :param partition_bin_count: Number of bins per partition
        :param partition_size: Size of each partition
        :param filters: Optional filters to apply during preparation
        :param obs_columns: Optional list of observation columns to include
        :param bucket_name: Optional GCS bucket name for metadata
        :param extract_bucket_path: Optional GCS path for metadata

        :raise google.api_core.exceptions.GoogleAPIError: If table preparation fails
        """
        prepare_extract_tables(
            client=self.client,
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=extract_table_prefix,
            features=features,
            extract_bin_size=extract_bin_size,
            assign_bin_by_category=assign_bin_by_category,
            extract_bin_category_column_name=extract_bin_category_column_name,
            random_seed_offset=random_seed_offset,
            partition_bin_count=partition_bin_count,
            partition_size=partition_size,
            filters=filters,
            obs_columns=obs_columns,
            bucket_name=bucket_name,
            extract_bucket_path=extract_bucket_path,
        )

    def extract_data(
        self,
        *,
        extract_table_prefix: str,
        start_bin: int,
        end_bin: int,
        output_dir: Path,
        obs_columns: list[str] | None = None,
        max_workers: int | None = None,
    ) -> None:
        """
        Extract data from prepared extract tables into AnnData files.

        :param extract_table_prefix: Prefix for extract table names
        :param start_bin: Starting bin number (inclusive)
        :param end_bin: Ending bin number (inclusive)
        :param output_dir: Local directory to save AnnData files
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
            start_bin=start_bin,
            end_bin=end_bin,
            output_dir=output_dir,
            obs_columns=obs_columns,
            max_workers=max_workers,
        )
