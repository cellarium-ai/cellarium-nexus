"""
Control and manage SOMA extract operations for Nexus.
"""

import logging
from typing import Any, Literal

from cellarium.nexus.clients import NexusBackendAPIClient
from cellarium.nexus.omics_datastore.soma_ops import TileDBSOMADataOperator
from cellarium.nexus.shared.schemas.omics_datastore import (
    GroupedCurriculumMetadata,
    RandomizedCurriculumMetadata,
)
from cellarium.nexus.shared.utils.workspace_file_manager import WorkspaceFileManager

logger = logging.getLogger(__name__)


class SomaDataOpsCoordinator:
    """
    Control and manage SOMA extract operations for Nexus.

    Provide a high-level interface for SOMA extract operations,
    delegating low-level operations to TileDBSOMADataOperator.
    """

    def __init__(self, *, experiment_uri: str, nexus_backend_api_url: str, bucket_name: str) -> None:
        """
        Initialize the SOMA data operations coordinator.

        :param experiment_uri: URI of the SOMA experiment
        :param nexus_backend_api_url: URL for Nexus backend API
        :param bucket_name: GCS bucket for storing curriculum metadata objects and extract outputs

        :raise ValueError: If experiment_uri or bucket_name is empty
        """
        if not experiment_uri:
            raise ValueError("experiment_uri cannot be empty")
        if not bucket_name:
            raise ValueError("bucket_name cannot be empty")

        self.experiment_uri = experiment_uri
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url)
        self.soma_operator = TileDBSOMADataOperator(experiment_uri=experiment_uri)
        self.workspace = WorkspaceFileManager(bucket_name=bucket_name)

        logger.info(f"Initialized SomaDataOpsCoordinator for {experiment_uri}")

    def _load_metadata_from_bucket(
        self,
        *,
        curriculum_metadata_path: str,
        extract_type: Literal["randomized", "grouped"] = "randomized",
    ) -> RandomizedCurriculumMetadata | GroupedCurriculumMetadata:
        """
        Load a SOMA extract metadata from cloud storage.

        :param curriculum_metadata_path: Path to the metadata JSON file within the bucket
        :param extract_type: Type of extraction - "randomized" or "grouped"

        :raise IOError: If file cannot be read
        :raise ValueError: If file is not valid JSON or metadata schema

        :return: Loaded RandomizedCurriculumMetadata or GroupedCurriculumMetadata
        """
        curriculum_metadata_json = self.workspace.load_json_from_bucket(remote_path=curriculum_metadata_path)
        if extract_type == "grouped":
            return GroupedCurriculumMetadata.model_validate(obj=curriculum_metadata_json)
        return RandomizedCurriculumMetadata.model_validate(obj=curriculum_metadata_json)

    def _update_curriculum_with_error(
        self,
        *,
        extract_name: str,
        error_message: str,
    ) -> None:
        """
        Update curriculum status to FAILED with error message.

        :param extract_name: Name of the curriculum/extract
        :param error_message: Error message to record
        """
        try:
            self.backend_client.update_curriculum(name=extract_name, status="FAILED")
            logger.error(f"Marked curriculum '{extract_name}' as FAILED: {error_message}")
        except Exception as e:
            logger.error(f"Failed to update curriculum status: {e}")

    def prepare_soma_extract(
        self,
        *,
        extract_name: str,
        creator_id: int,
        curriculum_metadata_path: str,
        filters: dict[str, Any] | None = None,
        var_filter_column: str | None = None,
        var_filter_values: list[str] | None = None,
        obs_columns: list[str] | None = None,
        var_columns: list[str] | None = None,
        x_layer: str = "X",
        # Randomized extraction params
        range_size: int | None = None,
        extract_bin_size: int | None = None,
        shuffle_ranges: bool = True,
        # Grouped extraction params
        extract_bin_keys: list[str] | None = None,
    ) -> RandomizedCurriculumMetadata | GroupedCurriculumMetadata:
        """
        Prepare a SOMA extract by computing the plan and registering the curriculum.

        Compute the extract plan, save it to cloud storage, and register the curriculum
        in the backend. Supports both randomized and grouped extraction modes.

        :param extract_name: Name for this extract/curriculum.
        :param creator_id: ID of the user creating the extract.
        :param curriculum_metadata_path: Path within bucket where the plan will be saved.
        :param filters: Optional filter conditions in Nexus format.
        :param var_filter_column: Name of the var column to filter features by.
        :param var_filter_values: List of values to match in the var filter column.
        :param obs_columns: List of obs columns to include in extraction.
        :param var_columns: List of var columns to include in extraction.
        :param x_layer: Name of the SOMA X layer to read counts from.
        :param range_size: Target number of cells per range (randomized extraction).
        :param extract_bin_size: Target cells per output bin.
        :param shuffle_ranges: Whether to shuffle the joinid ranges (randomized extraction).
        :param extract_bin_keys: List of obs column names to group by (grouped extraction).

        :raise SomaPlanningError: If plan computation fails
        :raise IOError: If cloud storage operations fail
        :raise HTTPError: If backend API calls fail
        :raise ValueError: If required parameters are missing

        :return: SOMA extract plan (RandomizedCurriculumMetadata or GroupedCurriculumMetadata)
        """
        extract_type = "grouped" if extract_bin_keys else "randomized"
        logger.info(f"Preparing SOMA {extract_type} extract '{extract_name}'")

        try:
            # Step 1: Compute extract plan
            plan = self.soma_operator.prepare_curriculum_metadata(
                filters=filters,
                range_size=range_size,
                extract_bin_size=extract_bin_size,
                shuffle_ranges=shuffle_ranges,
                var_filter_column=var_filter_column,
                var_filter_values=var_filter_values,
                obs_columns=obs_columns,
                var_columns=var_columns,
                x_layer=x_layer,
                extract_bin_keys=extract_bin_keys,
                bin_size=extract_bin_size,
            )
            logger.info(f"Computed {extract_type} extract plan: {plan.num_bins} bins, {plan.total_cells} total cells")

            # Step 2: Save plan to cloud storage
            full_plan_path = self.workspace.save_json_to_bucket(
                data=plan.model_dump(),
                remote_path=curriculum_metadata_path,
            )
            logger.info(f"Saved extract plan to {full_plan_path}")

            # Step 3: Register curriculum in backend
            self.backend_client.register_curriculum(
                name=extract_name,
                creator_id=creator_id,
                extract_bin_size=extract_bin_size,
                filters_json=filters,
            )
            logger.info(f"Registered curriculum '{extract_name}'")

            # Step 4: Update curriculum with metadata path and status
            self.backend_client.update_curriculum(
                name=extract_name,
                status="EXTRACTING",
                metadata_file_path=full_plan_path,
            )
            logger.info(f"Updated curriculum '{extract_name}' status to EXTRACTING")
            return plan

        except Exception as e:
            logger.error(f"Failed to prepare SOMA extract: {e}")
            self._update_curriculum_with_error(
                extract_name=extract_name,
                error_message=str(e),
            )
            raise

    def run_soma_extract(
        self,
        *,
        extract_name: str,
        curriculum_metadata_path: str,
        extract_bucket_path: str,
        extract_type: Literal["randomized", "grouped"] = "randomized",
        partition_index: int = 0,
        output_format: Literal["h5ad", "zarr"] = "h5ad",
        max_workers_extract: int | None = None,
        # Randomized-specific params
        max_ranges_per_partition: int | None = None,
        max_workers_shuffle: int | None = None,
        # Grouped-specific params
        max_bins_per_partition: int | None = None,
    ) -> None:
        """
        Run SOMA data extraction for specified ranges or bins.

        Load the extract curriculum_metadata from cloud storage and extract data to output files.
        Routes to randomized or grouped extraction based on extract_type.

        :param extract_name: Name of the extract/curriculum (for error reporting)
        :param curriculum_metadata_path: Path to the extract curriculum_metadata JSON within the bucket
        :param extract_bucket_path: Path within bucket for output files
        :param extract_type: Type of extraction - "randomized" or "grouped"
        :param partition_index: Zero-based index of the worker/partition
        :param output_format: Output format - "h5ad" or "zarr"
        :param max_workers_extract: Maximum parallel workers for extraction
        :param max_ranges_per_partition: Partition block size for randomized extraction
        :param max_workers_shuffle: Maximum parallel workers for shuffling (randomized only)
        :param max_bins_per_partition: Partition block size for grouped extraction

        :raise SomaExtractError: If extraction fails
        :raise IOError: If file operations fail
        """
        logger.info(f"Running SOMA {extract_type} extract for '{extract_name}'")

        try:
            # Step 1: Load curriculum_metadata from cloud storage
            curriculum_metadata = self._load_metadata_from_bucket(
                curriculum_metadata_path=curriculum_metadata_path,
                extract_type=extract_type,
            )
            logger.info(f"Loaded extract curriculum_metadata. Extracting partition index {partition_index}")

            # Step 2: Extract data using workspace for file management
            output_path = f"{extract_bucket_path}/extract_files"

            with self.workspace.temp_workspace() as paths:
                output_dir = paths["output"]

                if extract_type == "grouped":
                    # Grouped extraction - direct, no shuffle stage
                    self.soma_operator.extract_grouped(
                        curriculum_metadata=curriculum_metadata,
                        output_dir=output_dir,
                        partition_index=partition_index,
                        max_bins_per_partition=max_bins_per_partition,
                        output_format=output_format,
                        max_workers=max_workers_extract,
                    )
                else:
                    # Randomized extraction - two-stage with shuffle
                    temp_dir = paths["root"] / "temp_extract"
                    temp_dir.mkdir(exist_ok=True)

                    self.soma_operator.extract_randomized(
                        curriculum_metadata=curriculum_metadata,
                        output_dir=output_dir,
                        partition_index=partition_index,
                        max_ranges_per_partition=max_ranges_per_partition,
                        output_format=output_format,
                        temp_dir=temp_dir,
                        max_workers_extract=max_workers_extract,
                        max_workers_shuffle=max_workers_shuffle,
                        cleanup_temp=True,
                    )

                # Step 3: Upload output files to cloud storage
                self.workspace.upload_directory_to_bucket(
                    local_path=output_dir,
                    remote_path=output_path,
                )
                logger.info(f"Uploaded extract files to {self.workspace.bucket_name}/{output_path}")

        except Exception as e:
            logger.error(f"Failed to run SOMA extract: {e}")
            self._update_curriculum_with_error(
                extract_name=extract_name,
                error_message=str(e),
            )
            raise

    def mark_soma_curriculum_as_finished(
        self,
        *,
        extract_name: str,
        curriculum_metadata_path: str,
        extract_bucket_path: str,
        extract_type: Literal["randomized", "grouped"] = "randomized",
    ) -> None:
        """
        Mark a SOMA curriculum as finished and update with final metadata.

        Load the extract curriculum_metadata to get metrics and update the curriculum status
        to SUCCEEDED with cell count and extract info.

        :param extract_name: Name of the curriculum to update
        :param curriculum_metadata_path: Path to the extract curriculum_metadata JSON within the bucket
        :param extract_bucket_path: Path within bucket where extract files are stored
        :param extract_type: Type of extraction - "randomized" or "grouped"

        :raise IOError: If cloud storage operations fail
        :raise HTTPError: If backend API calls fail
        """
        logger.info(f"Marking SOMA curriculum '{extract_name}' as finished")

        try:
            # Step 1: Load curriculum_metadata to get metrics
            curriculum_metadata = self._load_metadata_from_bucket(
                curriculum_metadata_path=curriculum_metadata_path,
                extract_type=extract_type,
            )

            # Step 2: Update curriculum with final status and metrics
            extract_files_path = f"gs://{self.workspace.bucket_name}/{extract_bucket_path}/extract_files"
            metadata_file_path = f"gs://{self.workspace.bucket_name}/{curriculum_metadata_path}"

            self.backend_client.update_curriculum(
                name=extract_name,
                status="SUCCEEDED",
                cell_count=curriculum_metadata.total_cells,
                extract_bin_count=curriculum_metadata.num_bins,
                extract_files_path=extract_files_path,
                metadata_file_path=metadata_file_path,
            )
            logger.info(
                f"Marked curriculum '{extract_name}' as SUCCEEDED "
                f"(cells: {curriculum_metadata.total_cells}, bins: {curriculum_metadata.num_bins})"
            )

        except Exception as e:
            logger.error(f"Failed to mark curriculum as finished: {e}")
            self._update_curriculum_with_error(
                extract_name=extract_name,
                error_message=str(e),
            )
            raise
