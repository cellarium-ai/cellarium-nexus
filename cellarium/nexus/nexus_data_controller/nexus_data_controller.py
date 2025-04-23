"""
Control and manage Nexus data operations.
"""

import datetime
import logging
import os
import pathlib
import tempfile
from typing import Any, Sequence

import fastavro
import smart_open
from google.cloud import bigquery
from nexus.clients import NexusBackendAPIClient
from nexus.omics_datastore.bq_ops.bq_datastore_controller import BQDatastoreController
from nexus.omics_datastore.bq_ops.ingest.create_ingest_files import optimized_read_anndata

from cellarium.nexus.shared import schemas, utils

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class NexusDataController:
    """
    Control and manage Nexus data operations.

    This controller provides a high-level interface for data operations,
    delegating BigQuery-specific operations to BQDatastoreController.
    """

    def __init__(self, *, project_id: str, nexus_backend_api_url: str, bigquery_dataset: str) -> None:
        """
        Initialize the Nexus data controller.

        :param project_id: GCP project ID
        :param nexus_backend_api_url: URL for Nexus backend API
        :param bigquery_dataset: Name of BigQuery dataset

        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery client initialization fails
        """
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url)
        self.bq_client = bigquery.Client(project=project_id)
        self.project_id = project_id
        self.bq_controller = BQDatastoreController(
            client=self.bq_client,
            project=self.project_id,
            dataset=bigquery_dataset,
        )

    def __update_ingest_info_with_error(self, *, ingest_id: int, error_message: str) -> None:
        """
        Update ingest info with error message.

        :param ingest_id: ID of the ingest to update
        :param error_message: Error message to store

        :raise Exception: If update fails
        """
        self.backend_client.update_ingest_metadata_extra(
            ingest_id=ingest_id, new_metadata_extra={"error_occurred": error_message}
        )
        self.backend_client.update_ingest_status(ingest_id=ingest_id, new_status="FAILED")

    def create_ingest_files(
        self,
        *,
        input_file_path: str,
        tag: str | None,
        bigquery_dataset: str,
        bucket_name: str,
        bucket_stage_dir: str,
        column_mapping: dict[str, Any] | None = None,
    ) -> tuple[int, str]:
        """
        Create ingest files and prepare them for ingestion.

        :param input_file_path: Path to input file
        :param tag: Optional tag for the ingest
        :param bigquery_dataset: Name of BigQuery dataset
        :param bucket_name: Name of GCS bucket
        :param bucket_stage_dir: Staging directory in bucket
        :param column_mapping: Optional dictionary containing obs and var column mappings

        :raise Exception: If file creation or upload fails

        :return: Tuple of (ingest_id, nexus_uuid)
        """
        input_file_bucket_name, input_file_bucket_path = utils.gcp.get_bucket_name_and_file_path_from_gc_path(
            full_gs_path=input_file_path
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = pathlib.Path(temp_dir)
            logger.info(f"Created temporary directory `{temp_dir_path}`")

            local_input_data_dir = temp_dir_path / "input_files"
            local_output_dir = temp_dir_path / "output_files"

            os.makedirs(local_input_data_dir, exist_ok=True)
            os.makedirs(local_output_dir, exist_ok=True)

            # Create ingest info on backend
            ingest_info_api_struct = self.backend_client.create_ingest_file_info(bigquery_dataset=bigquery_dataset)

            # Download file from bucket
            logger.info("Downloading file from Bucket...")
            local_input_data_path = pathlib.Path(local_input_data_dir) / "adata.h5ad"
            utils.gcp.download_file_from_bucket(
                bucket_name=input_file_bucket_name,
                source_blob_name=input_file_bucket_path,
                destination_file_name=local_input_data_path,
            )

            logger.info("Reading the file...")
            try:
                adata = optimized_read_anndata(input_file_path=local_input_data_path)
                total_cells = adata.n_obs
                total_features = adata.n_vars

                adata.file.close()
                del adata

                cell_info_start_index, cell_info_end_index = self.backend_client.reserve_indexes_cell_info(
                    batch_size=total_cells
                )
                feature_info_start_index, feature_info_end_index = self.backend_client.reserve_indexes_feature_info(
                    batch_size=total_features
                )

                # Use instance bq_controller
                ingest_job_result = self.bq_controller.create_ingest_files(
                    adata_file_path=local_input_data_path,
                    tag=tag,
                    cell_info_start_index=cell_info_start_index,
                    cell_info_end_index=cell_info_end_index,
                    feature_info_start_index=feature_info_start_index,
                    feature_info_end_index=feature_info_end_index,
                    ingest_id=ingest_info_api_struct.id,
                    output_dir=local_output_dir,
                    column_mapping=column_mapping,
                )

                # Update ingest_info "metadata_extra" on the backend
                uns_json = ingest_job_result["adata_uns_clean"]
                self.backend_client.update_ingest_metadata_extra(
                    ingest_id=ingest_info_api_struct.id, new_metadata_extra=uns_json
                )

                # Upload to staging bucket
                ingest_stage_dir = f"{bucket_stage_dir}"
                utils.gcp.transfer_directory_to_bucket(
                    bucket_name=bucket_name, local_directory_path=local_output_dir, prefix=ingest_stage_dir
                )

                return ingest_info_api_struct.id, ingest_info_api_struct.nexus_uuid

            except Exception as e:
                self.__update_ingest_info_with_error(ingest_id=ingest_info_api_struct.id, error_message=str(e))
                raise

    def _read_gcs_avro_file(self, gcs_uri: str) -> list[dict]:
        """
        Read records from an Avro file in GCS.

        :param gcs_uri: Full GCS URI to the Avro file
        :return: List of records as dictionaries
        :raise IOError: If file cannot be read
        """
        with smart_open.open(gcs_uri, "rb") as f:
            avro_reader = fastavro.reader(f)
            return [record for record in avro_reader]

    def ingest_data_to_bigquery(
        self,
        *,
        bucket_name: str,
        bucket_stage_dir: str,
    ) -> None:
        """
        Ingest data from GCS into BigQuery tables.

        :param bucket_name: GCS bucket name containing the data
        :param bucket_stage_dir: Directory in the bucket containing staged files

        :raise google.api_core.exceptions.GoogleAPIError: If ingestion fails
        :raise Exception: If the ingest fails or status update fails
        """
        # Read ingest ID from avro file in GCS
        gcs_uri = f"gs://{bucket_name}/{bucket_stage_dir}/ingest-info.avro"
        try:
            records = self._read_gcs_avro_file(gcs_uri)
            if not records:
                raise ValueError("No records found in ingest-info.avro")
            if "id" not in records[0]:
                raise ValueError("No 'id' field found in ingest-info.avro record")
            ingest_id = records[0]["id"]
        except (IOError, ValueError) as e:
            logger.error(f"Error reading ingest ID from {gcs_uri}: {e}")
            raise

        try:
            self.bq_controller.ingest_data(
                gcs_bucket_name=bucket_name,
                gcs_stage_dir=bucket_stage_dir,
            )
            self.backend_client.ingest_from_avro(
                stage_dir=bucket_stage_dir,
                ingest_id=ingest_id,
            )
            self.backend_client.update_ingest_status(
                ingest_id=ingest_id,
                new_status="SUCCEEDED",
                ingest_finish_timestamp=datetime.datetime.now(),
            )
        except Exception as e:
            self.backend_client.update_ingest_status(
                ingest_id=ingest_id,
                new_status="FAILED",
                ingest_finish_timestamp=datetime.datetime.now(),
            )
            # Re-raise the exception to stop execution
            raise

    def create_bigquery_dataset(self, *, bigquery_dataset: str, location: str = "US") -> str:
        """
        Create a new BigQuery dataset and initialize required tables.

        :param bigquery_dataset: Name of the dataset to create
        :param location: Geographic location for the dataset

        :return: Name of the created dataset

        :raise google.api_core.exceptions.GoogleAPIError: If dataset creation fails
        """
        return self.bq_controller.initialize_bigquery_resources(location=location)

    def prepare_extract_tables(
        self,
        *,
        extract_table_prefix: str,
        features: Sequence[schemas.FeatureSchema],
        filters: dict[str, Any] | None = None,
        obs_columns: list[str] | None = None,
        extract_bin_size: int = 10000,
        assign_bin_by_category: bool = False,
        extract_bin_category_column_name: str | None = None,
        random_seed_offset: int = 0,
        partition_bin_count: int = 40000,
        partition_size: int = 10,
        bucket_name: str | None = None,
        extract_bucket_path: str | None = None,
    ) -> None:
        """
        Prepare extract tables for data extraction.

        :param extract_table_prefix: Prefix for extract table names
        :param features: Sequence of feature schema objects
        :param filters: Optional query filters to apply
        :param obs_columns: Optional list of observation columns to include
        :param extract_bin_size: Size of cell bins
        :param assign_bin_by_category: Whether to bin by category
        :param extract_bin_category_column_name: Column name for category binning
        :param random_seed_offset: Offset for randomization
        :param partition_bin_count: Number of partitions
        :param partition_size: Size of each partition
        :param bucket_name: GCS bucket name for metadata storage
        :param extract_bucket_path: Path within bucket for metadata storage

        :raise ValueError: If binning parameters are invalid
        :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
        :raise IOError: If metadata file operations fail
        """
        self.bq_controller.prepare_extract_tables(
            extract_table_prefix=extract_table_prefix,
            features=features,
            filters=filters,
            obs_columns=obs_columns,
            extract_bin_size=extract_bin_size,
            assign_bin_by_category=assign_bin_by_category,
            extract_bin_category_column_name=extract_bin_category_column_name,
            random_seed_offset=random_seed_offset,
            partition_bin_count=partition_bin_count,
            partition_size=partition_size,
            bucket_name=bucket_name,
            extract_bucket_path=extract_bucket_path,
        )

    def extract_data(
        self,
        *,
        extract_table_prefix: str,
        bins: list[int],
        bucket_name: str,
        extract_bucket_path: str,
        obs_columns: list[str] | None = None,
        max_workers: int | None = None,
    ) -> None:
        """
        Extract data from prepared extract tables into AnnData files and upload to GCS.

        :param extract_table_prefix: Prefix for extract table names
        :param bins: List of bin numbers to extract
        :param bucket_name: GCS bucket name
        :param extract_bucket_path: Path within bucket
        :param obs_columns: Optional list of observation columns to include
        :param max_workers: Maximum number of parallel workers

        :raise google.api_core.exceptions.GoogleAPIError: If extraction fails
        :raise IOError: If file operations fail
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = pathlib.Path(temp_dir)
            logger.info(f"Created temporary directory `{temp_dir_path}`")

            # Extract data locally
            self.bq_controller.extract_data(
                extract_table_prefix=extract_table_prefix,
                bins=bins,
                output_dir=temp_dir_path,
                obs_columns=obs_columns,
                max_workers=max_workers,
            )

            # Upload extracted files to GCS
            utils.gcp.transfer_directory_to_bucket(
                bucket_name=bucket_name,
                local_directory_path=temp_dir_path,
                prefix=f"{extract_bucket_path}/extract_files",
            )
