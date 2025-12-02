"""
Control and manage SOMA extract operations for Nexus.
"""

import logging
from typing import Any, Literal

from cellarium.nexus.clients import NexusBackendAPIClient
from cellarium.nexus.omics_datastore.soma_ops import TileDBSOMADataOperator
from cellarium.nexus.shared.schemas.omics_datastore import SomaExtractPlan
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
        :param bucket_name: GCS bucket for storing plans and extract outputs

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

    def _load_plan_from_bucket(self, *, plan_path: str) -> SomaExtractPlan:
        """
        Load a SOMA extract plan from cloud storage.

        :param plan_path: Path to the plan JSON file within the bucket

        :raise IOError: If file cannot be read
        :raise ValueError: If file is not valid JSON or plan schema

        :return: Loaded SomaExtractPlan
        """
        plan_data = self.workspace.load_json_from_bucket(remote_path=plan_path)
        return SomaExtractPlan.model_validate(plan_data)

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
        extract_bucket_path: str,
        plan_path: str,
        range_size: int,
        output_chunk_size: int,
        filters: dict[str, Any] | None = None,
        shuffle_ranges: bool = True,
        var_filter_column: str | None = None,
        var_filter_values: list[str] | None = None,
        obs_columns: list[str] | None = None,
        var_columns: list[str] | None = None,
        x_layer: str = "X",
    ) -> None:
        """
        Prepare a SOMA extract by computing the plan and registering the curriculum.

        Compute the extract plan (joinid ranges), save it to cloud storage, and register
        the curriculum in the backend.

        :param extract_name: Name for this extract/curriculum
        :param creator_id: ID of the user creating the extract
        :param extract_bucket_path: Path within bucket for this extract
        :param plan_path: Path within bucket where the plan will be saved
        :param range_size: Target number of cells per range
        :param output_chunk_size: Target cells per output chunk (for shuffling)
        :param filters: Optional filter conditions in Nexus format
        :param shuffle_ranges: Whether to shuffle the joinid ranges
        :param var_filter_column: Name of the var column to filter features by
        :param var_filter_values: List of values to match in the var filter column
        :param obs_columns: List of obs columns to include in extraction
        :param var_columns: List of var columns to include in extraction
        :param x_layer: Name of the SOMA X layer to read counts from

        :raise SomaPlanningError: If plan computation fails
        :raise IOError: If cloud storage operations fail
        :raise HTTPError: If backend API calls fail
        """
        logger.info(f"Preparing SOMA extract '{extract_name}' with range_size={range_size}")

        try:
            # Step 1: Compute extract plan
            plan = self.soma_operator.compute_extract_plan(
                filters=filters,
                range_size=range_size,
                output_chunk_size=output_chunk_size,
                shuffle_ranges=shuffle_ranges,
                var_filter_column=var_filter_column,
                var_filter_values=var_filter_values,
                obs_columns=obs_columns,
                var_columns=var_columns,
                x_layer=x_layer,
            )
            logger.info(f"Computed extract plan: {len(plan.joinid_ranges)} ranges, {plan.total_cells} total cells")

            # Step 2: Save plan to cloud storage
            full_plan_path = self.workspace.save_json_to_bucket(
                data=plan.model_dump(),
                remote_path=plan_path,
            )
            logger.info(f"Saved extract plan to {full_plan_path}")

            # Step 3: Register curriculum in backend
            self.backend_client.register_curriculum(
                name=extract_name,
                creator_id=creator_id,
                extract_bin_size=range_size,
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
        plan_path: str,
        extract_bucket_path: str,
        range_indices: list[int] | None = None,
        output_format: Literal["h5ad", "zarr"] = "h5ad",
        max_workers_extract: int | None = None,
        max_workers_shuffle: int | None = None,
    ) -> None:
        """
        Run SOMA data extraction for specified ranges.

        Load the extract plan from cloud storage, optionally subset to specific range indices,
        and extract data to output files. All data specification (obs_columns, var_columns,
        var_joinids, x_layer) is taken from the plan.

        :param extract_name: Name of the extract/curriculum (for error reporting)
        :param plan_path: Path to the extract plan JSON within the bucket
        :param extract_bucket_path: Path within bucket for output files
        :param range_indices: Indices of ranges to extract (None = all ranges)
        :param output_format: Output format - "h5ad" or "zarr"
        :param max_workers_extract: Maximum parallel workers for extraction
        :param max_workers_shuffle: Maximum parallel workers for shuffling

        :raise SomaExtractError: If extraction fails
        :raise IOError: If file operations fail
        """
        logger.info(f"Running SOMA extract for '{extract_name}'")

        try:
            # Step 1: Load plan from cloud storage
            plan = self._load_plan_from_bucket(plan_path=plan_path)
            logger.info(f"Loaded extract plan with {len(plan.joinid_ranges)} ranges")

            # Step 2: Subset ranges if indices provided
            if range_indices is not None:
                subset_ranges = [plan.joinid_ranges[i] for i in range_indices]
                plan = SomaExtractPlan(
                    experiment_uri=plan.experiment_uri,
                    value_filter=plan.value_filter,
                    joinid_ranges=subset_ranges,
                    total_cells=sum(r.end - r.start + 1 for r in subset_ranges),
                    range_size=plan.range_size,
                    output_chunk_size=plan.output_chunk_size,
                    filters=plan.filters,
                    var_joinids=plan.var_joinids,
                    var_filter_column=plan.var_filter_column,
                    var_filter_values=plan.var_filter_values,
                    obs_columns=plan.obs_columns,
                    var_columns=plan.var_columns,
                    x_layer=plan.x_layer,
                )
                logger.info(f"Subset to {len(subset_ranges)} ranges (indices {range_indices[0]}-{range_indices[-1]})")

            # Step 3: Extract data using workspace for file management
            output_path = f"{extract_bucket_path}/extract_files"

            with self.workspace.temp_workspace() as paths:
                output_dir = paths["output"]
                temp_dir = paths["root"] / "temp_extract"
                temp_dir.mkdir(exist_ok=True)

                # Run extraction with shuffle pipeline
                self.soma_operator.extract_ranges_shuffled(
                    plan=plan,
                    output_dir=output_dir,
                    output_format=output_format,
                    temp_dir=temp_dir,
                    max_workers_extract=max_workers_extract,
                    max_workers_shuffle=max_workers_shuffle,
                    cleanup_temp=True,
                )

                # Step 4: Upload output files to cloud storage
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
        plan_path: str,
        extract_bucket_path: str,
    ) -> None:
        """
        Mark a SOMA curriculum as finished and update with final metadata.

        Load the extract plan to get metrics and update the curriculum status
        to SUCCEEDED with cell count and extract info.

        :param extract_name: Name of the curriculum to update
        :param plan_path: Path to the extract plan JSON within the bucket
        :param extract_bucket_path: Path within bucket where extract files are stored

        :raise IOError: If cloud storage operations fail
        :raise HTTPError: If backend API calls fail
        """
        logger.info(f"Marking SOMA curriculum '{extract_name}' as finished")

        try:
            # Step 1: Load plan to get metrics
            plan = self._load_plan_from_bucket(plan_path=plan_path)

            # Step 2: Update curriculum with final status and metrics
            extract_files_path = f"gs://{self.workspace.bucket_name}/{extract_bucket_path}/extract_files"
            metadata_file_path = f"gs://{self.workspace.bucket_name}/{plan_path}"

            self.backend_client.update_curriculum(
                name=extract_name,
                status="SUCCEEDED",
                cell_count=plan.total_cells,
                extract_bin_count=len(plan.joinid_ranges),
                extract_files_path=extract_files_path,
                metadata_file_path=metadata_file_path,
            )
            logger.info(
                f"Marked curriculum '{extract_name}' as SUCCEEDED "
                f"(cells: {plan.total_cells}, ranges: {len(plan.joinid_ranges)})"
            )

        except Exception as e:
            logger.error(f"Failed to mark curriculum as finished: {e}")
            self._update_curriculum_with_error(
                extract_name=extract_name,
                error_message=str(e),
            )
            raise
