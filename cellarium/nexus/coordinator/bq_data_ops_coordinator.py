"""
Control and manage BigQuery data operations for Nexus.
"""

import concurrent.futures
import datetime
import json
import logging
import os
import pathlib
import tempfile
from typing import Any, ContextManager, Sequence

import fastavro
import smart_open
from google.cloud import bigquery

from cellarium.nexus.clients import NexusBackendAPIClient
from cellarium.nexus.coordinator import constants, exceptions
from cellarium.nexus.omics_datastore.bq_ops import BigQueryDataOperator, BigQueryDataValidator
from cellarium.nexus.omics_datastore.bq_ops.ingest.create_ingest_files import optimized_read_anndata
from cellarium.nexus.shared import schemas, utils
from cellarium.nexus.shared.schemas.omics_datastore import ExtractMetadata

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class BQDataOpsCoordinator:
    """
    Control and manage BigQuery data operations for Nexus.

    Provide a high-level interface for BigQuery ingest and extract operations,
    delegating low-level operations to BigQueryDataOperator.
    """

    def __init__(self, *, project_id: str, nexus_backend_api_url: str, bigquery_dataset: str) -> None:
        """
        Initialize the BigQuery data operations coordinator.

        :param project_id: GCP project ID
        :param nexus_backend_api_url: URL for Nexus backend API
        :param bigquery_dataset: Name of BigQuery dataset

        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery client initialization fails
        """
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url)
        self.bq_client = bigquery.Client(project=project_id)
        self.project_id = project_id
        self.bq_data_operator = BigQueryDataOperator(
            client=self.bq_client,
            project=self.project_id,
            dataset=bigquery_dataset,
        )

    def _update_ingest_info_with_error(self, *, ingest_id: int, error_message: str | dict[str, str]) -> None:
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

    def _update_ingest_info_with_success(self, *, ingest_id: int) -> None:
        self.backend_client.update_ingest_status(
            ingest_id=ingest_id, new_status="SUCCEEDED", ingest_finish_timestamp=datetime.datetime.now()
        )

    def _run_validation_methods(
        self, *, local_input_data_path: pathlib.Path, validation_methods: list[str] | None, ingest_id: int
    ) -> None:
        if validation_methods:
            logger.info("Validating the file...")
            validation_results = BigQueryDataValidator.call_validation_methods(
                validation_methods=validation_methods, anndata_file_local_path=local_input_data_path
            )
            validation_messages = {}
            error_occurred = False

            for validation_method, (is_valid, messages, has_warning) in zip(validation_methods, validation_results):
                if not is_valid:
                    error_occurred = True
                    validation_messages[validation_method] = messages

            if error_occurred:
                logger.error(f"Validation messages: {validation_messages}")
                self._update_ingest_info_with_error(ingest_id=ingest_id, error_message=validation_messages)
                raise exceptions.NexusDataOpsValidationError("Validation failed for file")

    def create_ingest_files(
        self,
        *,
        input_file_path: str,
        tag: str | None,
        bigquery_dataset: str,
        bucket_name: str,
        bucket_stage_dir: str,
        max_input_data_size: int = constants.MAX_INPUT_INGEST_FILE_SIZE_BYTES,
        column_mapping: dict[str, Any] | None = None,
        uns_keys_to_keep: list[str] | None = None,
        validation_methods: list[str] | None = None,
        ingest_id: int | None = None,
        ingest_file_id: int | None = None,
    ) -> tuple[int, str]:
        """
        Create ingest files and prepare them for ingestion.

        :param input_file_path: Path to input file
        :param tag: Optional tag for the ingest
        :param bigquery_dataset: Name of BigQuery dataset
        :param bucket_name: Name of GCS bucket
        :param bucket_stage_dir: Staging directory in bucket
        :param max_input_data_size: Maximum size of input file in bytes. If the file is larger than this, an exception
            will be raised.
        :param column_mapping: Optional dictionary containing obs and var column mappings
        :param uns_keys_to_keep: Optional list of keys to keep in the `uns` JSON blob. If not provided,
            all keys will be kept
        :param validation_methods: Optional list of validation methods to run on the data. If not provided, none applied
        :param ingest_id: Optional parent ingest ID to associate with this file
        :param ingest_file_id: Optional ingest file ID to reuse

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

            # Create or reuse ingest info on backend
            if ingest_file_id is not None:
                ingest_info_api_struct = self.backend_client.update_ingest_status(
                    ingest_id=ingest_file_id, new_status="STARTED"
                )
            else:
                ingest_info_api_struct = self.backend_client.create_ingest_file_info(
                    omics_dataset=bigquery_dataset,
                    ingest_id=ingest_id,
                    gcs_file_path=input_file_path,
                    tag=tag,
                )

            logger.info("Validating the file size...")
            BigQueryDataValidator.validate_remote_file_size(
                adata_gcs_path=input_file_path, max_size_bytes=max_input_data_size
            )

            logger.info("Downloading file from Bucket...")
            local_input_data_path = pathlib.Path(local_input_data_dir) / "adata.h5ad"
            utils.gcp.download_file_from_bucket(
                bucket_name=input_file_bucket_name,
                source_blob_name=input_file_bucket_path,
                destination_file_name=local_input_data_path,
            )

            self._run_validation_methods(
                local_input_data_path=local_input_data_path,
                validation_methods=validation_methods,
                ingest_id=ingest_info_api_struct.id,
            )

            logger.info("Reading the file...")
            try:
                adata = optimized_read_anndata(input_file_path=local_input_data_path)
                total_cells = adata.n_obs
                total_features = adata.n_vars

                adata.file.close()
                del adata

                cell_info_start_index, cell_info_end_index = self.backend_client.reserve_indexes_cell_info(
                    bigquery_dataset=bigquery_dataset,
                    batch_size=total_cells,
                )
                feature_info_start_index, feature_info_end_index = self.backend_client.reserve_indexes_feature_info(
                    bigquery_dataset=bigquery_dataset,
                    batch_size=total_features,
                )

                ingest_job_result = self.bq_data_operator.create_ingest_files(
                    adata_file_path=local_input_data_path,
                    tag=tag,
                    cell_info_start_index=cell_info_start_index,
                    cell_info_end_index=cell_info_end_index,
                    feature_info_start_index=feature_info_start_index,
                    feature_info_end_index=feature_info_end_index,
                    ingest_id=ingest_info_api_struct.id,
                    output_dir=local_output_dir,
                    column_mapping=column_mapping,
                    uns_keys_to_keep=uns_keys_to_keep,
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
                self._update_ingest_info_with_error(ingest_id=ingest_info_api_struct.id, error_message=str(e))
                raise

    def _read_gcs_avro_file(self, gcs_uri: str) -> list[dict]:
        """
        Read an Avro file from GCS.

        :param gcs_uri: GCS URI of the Avro file to read

        :raise: IOError: If the file cannot be read
        :raise: ValueError: If the file is not a valid Avro file

        :return: List of records from the Avro file
        """
        with smart_open.open(gcs_uri, "rb") as f:
            avro_reader = fastavro.reader(f)
            return [record for record in avro_reader]

    def _ingest_data_to_bigquery(self, bucket_name: str, bucket_stage_dir: str) -> ContextManager:
        logger.info("Starting BigQuery data ingest...")
        return self.bq_data_operator.ingest_data(
            gcs_bucket_name=bucket_name,
            gcs_stage_dir=bucket_stage_dir,
        )

    @staticmethod
    def _cleanup_ingest_files(bucket_name: str, bucket_stage_dir: str) -> list[str]:
        """
        Clean up ingest files from a GCS bucket after successful ingestion.

        :param bucket_name: GCS bucket name containing the data
        :param bucket_stage_dir: Directory in the bucket containing staged files

        :raise: google.cloud.exceptions.GoogleAPIError: If deletion fails

        :return: List of deleted blob names
        """
        logger.info(f"Cleaning up ingest files from gs://{bucket_name}/{bucket_stage_dir}")
        deleted_files = utils.gcp.delete_files_from_bucket(bucket_name=bucket_name, prefix=bucket_stage_dir)
        logger.info(f"Deleted {len(deleted_files)} files from gs://{bucket_name}/{bucket_stage_dir}")
        return deleted_files

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
            with self._ingest_data_to_bigquery(bucket_name=bucket_name, bucket_stage_dir=bucket_stage_dir):
                # Nothing to do inside; commit happens on clean exit of this block.
                pass

            self._cleanup_ingest_files(bucket_name=bucket_name, bucket_stage_dir=bucket_stage_dir)
        except Exception as e:
            self._update_ingest_info_with_error(ingest_id=ingest_id, error_message=str(e))
            # Re-raise the exception to stop execution
            raise
        else:
            self._update_ingest_info_with_success(ingest_id=ingest_id)

    def ingest_data_to_bigquery_parallel(
        self,
        *,
        bucket_name: str,
        bucket_stage_dirs: list[str],
        num_workers: int = 4,
    ) -> None:
        """
        Ingest data from multiple GCS stage directories into BigQuery tables asynchronously.

        :param bucket_name: GCS bucket name containing the data
        :param bucket_stage_dirs: List of directories in the bucket containing staged files
        :param num_workers: Number of concurrent workers for processing

        :raise: ValueError: If bucket_stage_dirs is empty
        :raise: Exception: If any of the ingestion operations fail (logged but not propagated)
        """
        if not bucket_stage_dirs:
            raise ValueError("bucket_stage_dirs cannot be empty")

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            # Submit all tasks to the executor
            future_to_dir = {}
            for stage_dir in bucket_stage_dirs:
                future = executor.submit(
                    self.ingest_data_to_bigquery,
                    bucket_name=bucket_name,
                    bucket_stage_dir=stage_dir,
                )
                future_to_dir[future] = stage_dir

            # Process results as they complete
            error_happened = False
            for future in concurrent.futures.as_completed(future_to_dir):
                stage_dir = future_to_dir[future]
                try:
                    future.result()
                    logger.info(f"Successfully ingested data from {stage_dir}")
                except Exception as exc:
                    error_happened = True
                    logger.error(f"Failed to ingest data from {stage_dir}: {exc}")

        logger.info("Resetting Backend cache...")
        backend_reset_result = self.backend_client.reset_backend_cache()
        logger.info(f"Reset cache keys: {backend_reset_result.repopulated_keys}")

        if error_happened:
            raise exceptions.NexusDataOpsIngestError("Some ingest operations failed")

    def prepare_extract_tables(
        self,
        *,
        extract_name: str,
        features: Sequence[schemas.FeatureSchema],
        categorical_column_count_limit: int,
        creator_id: int,
        bucket_name: str,
        extract_bucket_path: str,
        extract_bin_size: int = 10000,
        random_seed_offset: int = 0,
        partition_bin_count: int = 40000,
        partition_size: int = 10,
        extract_bin_keys: list[str] | None = None,
        filters: dict[str, Any] | None = None,
        obs_columns: list[str] | None = None,
        metadata_extra_columns: list[str] | None = None,
    ) -> None:
        """
        Prepare extract tables for data extraction.

        :param extract_name: Prefix for extract table names
        :param features: Sequence of feature schema objects
        :param categorical_column_count_limit: Maximum number of categories per categorical column to be considered as
            categorical. If the number of categories exceeds this limit, the column will not be unified across all
            extract files.
        :param creator_id: ID of the curriculum creator
        :param bucket_name: GCS bucket name for metadata storage
        :param extract_bucket_path: Path within bucket for metadata storage
        :param extract_bin_size: Size of cell bins
        :param random_seed_offset: Offset for randomization
        :param partition_bin_count: Number of partitions
        :param partition_size: Size of each partition
        :param extract_bin_keys: Optional list of keys to use for binning. If not provided, all keys will be used.
        :param filters: Optional query filters to apply
        :param obs_columns: Optional list of observation columns to include
        :param metadata_extra_columns: Optional list of metadata extra columns to include from `metadata_extra` JSON
            blob. If not provided, none will be included.

        :raise ValueError: If binning parameters are invalid
        :raise google.api_core.exceptions.GoogleAPIError: If table creation fails
        :raise IOError: If metadata file operations fail
        """
        backend_client = NexusBackendAPIClient(api_url=self.backend_client.api_url)
        backend_client.register_curriculum(
            name=extract_name, creator_id=creator_id, extract_bin_size=extract_bin_size, filters_json=filters
        )
        try:
            extract_metadata = self.bq_data_operator.prepare_extract_tables(
                extract_table_prefix=extract_name,
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

            # Save metadata to GCS
            metadata_path = f"{extract_bucket_path}/{constants.EXTRACT_METADATA_FILE_NAME}"

            # Create a temporary local file
            with tempfile.NamedTemporaryFile(mode="w", suffix=".json") as temp_file:
                json.dump(extract_metadata.model_dump(), temp_file, indent=2)
                temp_file.flush()

                # Upload to GCS
                logger.info("Uploading metadata to GCS")
                utils.gcp.upload_file_to_bucket(
                    local_file_path=temp_file.name, bucket_name=bucket_name, blob_name=metadata_path
                )

            # Update curriculum with metadata file path
            backend_client.update_curriculum(name=extract_name, metadata_file_path=metadata_path, status="EXTRACTING")
        except Exception:
            backend_client.update_curriculum(name=extract_name, status="FAILED")
            raise

    def extract_data(
        self,
        *,
        extract_name: str,
        bins: list[int],
        bucket_name: str,
        extract_bucket_path: str,
        obs_columns: list[str] | None = None,
        max_workers: int | None = None,
    ) -> None:
        """
        Extract data from prepared extract tables into AnnData files and upload to GCS.

        :param extract_name: Prefix for extract table names
        :param bins: List of bin numbers to extract
        :param bucket_name: GCS bucket name
        :param extract_bucket_path: Path within bucket
        :param obs_columns: Optional list of observation columns to include
        :param max_workers: Maximum number of parallel workers

        :raise google.api_core.exceptions.GoogleAPIError: If extraction fails
        :raise IOError: If file operations fail
        :raise ValueError: If extract metadata is invalid
        """
        # Load extract metadata from GCS
        metadata_path = f"gs://{bucket_name}/{extract_bucket_path}/{constants.EXTRACT_METADATA_FILE_NAME}"
        with smart_open.open(metadata_path, "r") as f:
            metadata_dict = json.load(f)
            extract_metadata = schemas.ExtractMetadata(**metadata_dict)
        logger.info(f"Loaded extract metadata from `{metadata_path}`")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = pathlib.Path(temp_dir)
            logger.info(f"Created temporary directory `{temp_dir_path}`")

            # Extract data locally
            self.bq_data_operator.extract_data(
                extract_table_prefix=extract_name,
                bins=bins,
                output_dir=temp_dir_path,
                extract_metadata=extract_metadata,
                obs_columns=obs_columns,
                max_workers=max_workers,
            )

            # Upload extracted files to GCS
            utils.gcp.transfer_directory_to_bucket(
                bucket_name=bucket_name,
                local_directory_path=temp_dir_path,
                prefix=f"{extract_bucket_path}/extract_files",
            )

    def mark_curriculum_as_finished(
        self,
        *,
        extract_name: str,
        bucket_name: str,
        extract_bucket_path: str,
    ) -> None:
        """
        Mark a curriculum as finished and succeeded, including extract metadata.

        Retrieves metadata information from the metadata file and updates the curriculum
        with this information, including the number of bins, extract files path, and
        metadata file path. If any error occurs during the process, the curriculum is marked as failed.

        :param extract_name: Name of the curriculum to mark as finished
        :param bucket_name: GCS bucket name containing the extract files
        :param extract_bucket_path: Path within bucket for the extract

        :raise Exception: If any error occurs during the process (after marking curriculum as failed)
        """
        # Define paths
        extract_files_path = f"{extract_bucket_path}/extract_files"
        metadata_path = f"{extract_bucket_path}/{constants.EXTRACT_METADATA_FILE_NAME}"
        gcs_metadata_path = f"gs://{bucket_name}/{metadata_path}"

        try:
            # Read metadata from the JSON file
            logger.info(f"Reading metadata from {gcs_metadata_path}")
            with smart_open.open(gcs_metadata_path, "r") as f:
                metadata_dict = json.load(f)

            # Parse metadata using the Pydantic model
            extract_metadata = ExtractMetadata.model_validate(metadata_dict)

            # Get required information
            total_bins = extract_metadata.total_bins

            # Calculate cell count using the model's method
            cell_count = extract_metadata.calculate_cell_count()

            logger.info(f"Curriculum {extract_name} has {total_bins} bins and {cell_count} cells")

            # Update curriculum with metadata and mark as succeeded
            self.backend_client.update_curriculum(
                name=extract_name,
                cell_count=cell_count,
                extract_bin_count=total_bins,
                extract_files_path=extract_files_path,
                metadata_file_path=metadata_path,
                status="SUCCEEDED",
            )

            logger.info(f"Successfully marked curriculum {extract_name} as finished")

        except Exception as e:
            error_message = f"Failed to mark curriculum {extract_name} as finished: {e}"
            logger.error(error_message)

            # Mark the curriculum as failed
            try:
                self.backend_client.update_curriculum(name=extract_name, status="FAILED")
                logger.info(f"Marked curriculum {extract_name} as failed")
            except Exception as api_error:
                logger.error(f"Failed to mark curriculum as failed: {api_error}")

            # Re-raise the original exception
            raise
