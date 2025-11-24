"""
TileDB SOMA data operator for Nexus omics datastore.

This module provides the main operator class for interacting with SOMA experiments.
"""

import logging
import shutil
import tempfile
from pathlib import Path
from typing import Literal

import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaReadError
from cellarium.nexus.omics_datastore.soma_ops.extract import extract_ranges, shuffle_extracted_chunks
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter
from cellarium.nexus.omics_datastore.soma_ops.planning import plan_soma_extract
from cellarium.nexus.shared.schemas.omics_datastore import SomaExtractPlan

logger = logging.getLogger(__name__)


class TileDBSOMADataOperator:
    """
    Control and manage TileDB SOMA operations for the Nexus datastore.

    Provide low-level operations for interacting with a SOMA experiment:
    - Count cells with filters
    - Compute extract plans (soma_joinid ranges)
    - Extract data for joinid ranges into AnnData
    """

    def __init__(
        self,
        *,
        experiment_uri: str,
    ) -> None:
        """
        Initialize the SOMA datastore operator.

        :param experiment_uri: URI of the SOMA experiment

        :raise ValueError: If experiment_uri is invalid
        """
        if not experiment_uri:
            raise ValueError("experiment_uri cannot be empty")

        self.experiment_uri = experiment_uri

        logger.info(f"Initialized TileDBSOMADataOperator for {experiment_uri}")

    def count_cells(
        self,
        *,
        filters: dict[str, object] | None = None,
    ) -> int:
        """
        Count cells in the SOMA obs table matching the given filters.

        :param filters: Dict of filter conditions using the Nexus format

        :raise SomaReadError: If SOMA read fails

        :return: Total number of matching cells
        """
        try:
            value_filter = build_soma_value_filter(filters=filters)

            logger.info(f"Counting cells with filter: {value_filter if value_filter else 'no filter'}")

            with tiledbsoma.open(self.experiment_uri, mode="r") as exp:
                # Read only soma_joinid column with the filter
                obs_query = exp.obs.read(
                    column_names=["soma_joinid"],
                    value_filter=value_filter if value_filter else None,
                )

                # Get the count
                obs_df = obs_query.concat().to_pandas()
                count = len(obs_df)

                logger.info(f"Found {count} cells matching the filter")
                return count

        except Exception as e:
            logger.error(f"Failed to count cells: {e}")
            raise SomaReadError("SOMA cell count operation failed") from e

    def compute_extract_plan(
        self,
        *,
        filters: dict[str, object] | None,
        range_size: int,
        output_chunk_size: int | None = None,
        shuffle_ranges: bool = True,
    ) -> SomaExtractPlan:
        """
        Compute a SOMA extract plan from filter dict.

        Delegate to the planning module to compute contiguous joinid ranges.

        :param filters: Dict of filter conditions using the Nexus format
        :param range_size: Target number of cells per range (for extraction)
        :param output_chunk_size: Target number of cells per output chunk (for shuffling, default: range_size)
        :param shuffle_ranges: Whether to shuffle the resulting joinid ranges

        :raise SomaPlanningError: If SOMA reads fail
        :raise SomaReadError: If SOMA reads fail
        :raise ValueError: If range_size is not positive

        :return: SomaExtractPlan with joinid ranges and total cell count
        """
        return plan_soma_extract(
            experiment_uri=self.experiment_uri,
            filters=filters,
            range_size=range_size,
            output_chunk_size=output_chunk_size,
            shuffle_ranges=shuffle_ranges,
        )

    def extract_ranges_to_anndata(
        self,
        *,
        plan: SomaExtractPlan,
        output_dir: Path,
        obs_columns: list[str] | None = None,
        var_columns: list[str] | None = None,
        x_layer: str = "X",
        output_format: Literal["zarr", "h5ad"] = "h5ad",
        max_workers: int | None = None,
    ) -> None:
        """
        Extract SOMA data for the given joinid ranges into AnnData files.

        Delegate to the extract module to perform parallel extraction.

        :param plan: SOMA extract plan with joinid ranges and value_filter
        :param output_dir: Local directory to save AnnData files
        :param obs_columns: Optional list of obs columns to include
        :param var_columns: Optional list of var columns to include
        :param x_layer: Name of the SOMA X layer to read counts from
        :param output_format: Output format - "zarr" or "h5ad" (default: "h5ad")
        :param max_workers: Maximum number of parallel workers

        :raise SomaExtractError: If SOMA reads fail
        :raise IOError: If file operations fail
        """
        logger.info(f"Extracting {len(plan.joinid_ranges)} ranges to {output_dir} (format: {output_format})")

        extract_ranges(
            plan=plan,
            output_dir=output_dir,
            obs_columns=obs_columns,
            var_columns=var_columns,
            x_layer=x_layer,
            output_format=output_format,
            max_workers=max_workers,
        )

        logger.info("Extract operation completed successfully")

    def extract_ranges_shuffled(
        self,
        *,
        plan: SomaExtractPlan,
        output_dir: Path,
        output_chunk_size: int | None = None,
        output_format: Literal["zarr", "h5ad"] = "h5ad",
        temp_dir: Path | None = None,
        obs_columns: list[str] | None = None,
        var_columns: list[str] | None = None,
        x_layer: str = "X",
        max_workers_extract: int | None = None,
        max_workers_shuffle: int | None = None,
        cleanup_temp: bool = True,
    ) -> None:
        """
        Extract and shuffle cells across chunks.

        Two-stage process for true cell randomization:
        1. Extract contiguous ranges to temp directory in Zarr format (fast SOMA reads)
        2. Shuffle cells across final output chunks (true randomization)

        Delegate to extract and shuffle modules for implementation.

        :param plan: SOMA extract plan with joinid ranges and value_filter
        :param output_dir: Final output directory for shuffled chunks
        :param output_chunk_size: Number of cells per output chunk (default: plan.range_size)
        :param output_format: Output format - "zarr" or "h5ad" (default: "h5ad")
        :param temp_dir: Temporary directory for contiguous extracts (auto-created if None)
        :param obs_columns: Optional list of obs columns to include
        :param var_columns: Optional list of var columns to include
        :param x_layer: Name of the SOMA X layer to read counts from
        :param max_workers_extract: Maximum parallel workers for extraction (network I/O intensive)
        :param max_workers_shuffle: Maximum parallel workers for shuffling (CPU/memory intensive)
        :param cleanup_temp: Whether to delete temp directory after shuffling

        :raise SomaExtractError: If SOMA reads fail
        :raise IOError: If file operations fail
        :raise ValueError: If output_format is invalid
        """
        # Use plan.output_chunk_size if set, otherwise fall back to output_chunk_size param or plan.range_size
        final_chunk_size = plan.output_chunk_size or output_chunk_size or plan.range_size

        logger.info(
            f"Extracting and shuffling {len(plan.joinid_ranges)} ranges to {output_dir} "
            f"(output_chunk_size: {final_chunk_size}, format: {output_format})"
        )

        # Stage 1: Extract contiguous ranges to temp directory
        if temp_dir is None:
            temp_dir = Path(tempfile.mkdtemp(prefix="soma_extract_temp_"))
            logger.info(f"Created temporary directory: {temp_dir}")
        else:
            temp_dir.mkdir(parents=True, exist_ok=True)

        try:
            logger.info("Stage 1: Extracting contiguous ranges to Zarr format...")
            extract_ranges(
                plan=plan,
                output_dir=temp_dir,
                obs_columns=obs_columns,
                var_columns=var_columns,
                x_layer=x_layer,
                output_format="zarr",  # Use Zarr for temp files (fast random access)
                max_workers=max_workers_extract,
            )

            # Stage 2: Shuffle cells across chunks
            logger.info("Stage 2: Shuffling cells across chunks...")
            shuffle_extracted_chunks(
                input_dir=temp_dir,
                output_dir=output_dir,
                chunk_size=final_chunk_size,
                input_format="zarr",  # Read from Zarr temp files
                output_format=output_format,  # Write in requested format
                max_workers=max_workers_shuffle,
            )

            logger.info("Extract and shuffle operation completed successfully")

        finally:
            # Cleanup temp files
            if cleanup_temp and temp_dir.exists():
                logger.info(f"Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
