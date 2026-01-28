"""
TileDB SOMA data operator for Nexus omics datastore.

This module provides the main operator class for interacting with SOMA experiments.
"""

import logging
from pathlib import Path
from typing import Any, Literal

from cellarium.nexus.omics_datastore.protocols import DataOperatorProtocol
from cellarium.nexus.omics_datastore.soma_ops import _extract, _queries
from cellarium.nexus.shared.schemas.omics_datastore import (
    GroupedCurriculumMetadata,
    RandomizedCurriculumMetadata,
)

logger = logging.getLogger(__name__)


class TileDBSOMADataOperator(DataOperatorProtocol):
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

    def count_cells(self, *, filter_statements: dict[str, Any] | None = None) -> int:
        """
        Count cells in the SOMA obs table matching the given filters.

        :param filter_statements: Dict of filter conditions using the Nexus format

        :raise SomaReadError: If SOMA read fails

        :return: Total number of matching cells
        """
        return _queries.count_obs(
            experiment_uri=self.experiment_uri,
            filter_statements=filter_statements,
        )

    def get_categorical_obs_columns(self, *, threshold: int, exclude: list[str] | None = None) -> set[str]:
        """
        Determine categorical string columns in obs by distinct count threshold.

        :param threshold: Maximum distinct values to consider a column categorical
        :param exclude: Optional list of column names to skip

        :raise SomaReadError: If SOMA read fails

        :return: Set of column names that are categorical
        """
        return _queries.get_obs_string_columns(
            experiment_uri=self.experiment_uri,
            exclude=exclude,
        )

    def get_distinct_obs_values(self, *, column_name: str) -> list[str]:
        """
        Fetch distinct values for an obs column.

        :param column_name: Column to retrieve values for

        :raise SomaReadError: If SOMA read fails

        :return: List of distinct values
        """
        return _queries.get_obs_distinct_values(
            experiment_uri=self.experiment_uri,
            column_name=column_name,
        )

    def prepare_curriculum_metadata(
        self,
        *,
        filters: dict[str, object] | None,
        range_size: int | None = None,
        extract_bin_size: int | None = None,
        shuffle_ranges: bool = True,
        var_filter_column: str | None = None,
        var_filter_values: list[str] | None = None,
        obs_columns: list[str] | None = None,
        var_columns: list[str] | None = None,
        x_layer: str = "X",
        extract_bin_keys: list[str] | None = None,
        bin_size: int | None = None,
    ) -> RandomizedCurriculumMetadata | GroupedCurriculumMetadata:
        """
        Compute a SOMA extract plan from filter dict.

        If extract_bin_keys is provided, compute grouped bins where cells from the same
        group stay together. Otherwise, compute contiguous joinid ranges for randomized
        extraction.

        :param filters: Dict of filter conditions using the Nexus format.
        :param range_size: Target number of cells per range (for randomized extraction).
        :param extract_bin_size: Target number of cells per extract bin (for randomized shuffling).
        :param shuffle_ranges: Whether to shuffle the resulting joinid ranges (randomized only).
        :param var_filter_column: Name of the var column to filter features by.
        :param var_filter_values: List of values to match in the var filter column.
        :param obs_columns: List of obs columns to include in extraction.
        :param var_columns: List of var columns to include in extraction.
        :param x_layer: Name of the SOMA X layer to read counts from.
        :param extract_bin_keys: List of obs column names to group by (for grouped extraction).
        :param bin_size: Maximum number of cells per bin (for grouped extraction).

        :raise SomaPlanningError: If SOMA reads fail
        :raise SomaReadError: If SOMA reads fail
        :raise ValueError: If required parameters are missing or invalid

        :return: RandomizedCurriculumMetadata or GroupedCurriculumMetadata
        """
        return _extract.plan_curriculum(
            experiment_uri=self.experiment_uri,
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
            bin_size=bin_size,
        )

    def extract_randomized(
        self,
        *,
        curriculum_metadata: RandomizedCurriculumMetadata,
        output_dir: Path,
        partition_index: int = 0,
        max_ranges_per_partition: int | None = None,
        output_format: Literal["zarr", "h5ad"] = "h5ad",
        temp_dir: Path | None = None,
        max_workers_extract: int | None = None,
        max_workers_shuffle: int | None = None,
        cleanup_temp: bool = True,
        verbose: bool = False,
    ) -> None:
        """
        Extract and shuffle cells across chunks.

        Two-stage process for true cell randomization:
        1. Extract contiguous ranges to temp directory in Zarr format (fast SOMA reads)
        2. Shuffle cells across final output chunks (true randomization)

        Delegate to _extract and shuffle modules for implementation.
        All data specification (obs_columns, var_columns, var_joinids, x_layer)
        is taken from the plan.

        :param curriculum_metadata: SOMA curriculum metadata with all data specification.
        :param output_dir: Final output directory for shuffled chunks.
        :param partition_index: Index used for slicing ranges and output chunk indexes.
            Needed for distributing extracting over multiple distributed VMs. Default is 0 (single VM execution).
        :param max_ranges_per_partition: Partition block size. Default is None, this means it will use all ranges
        :param output_format: Output format - "zarr" or "h5ad" (default: "h5ad").
        :param temp_dir: Temporary directory for contiguous extracts (auto-created if None).
        :param max_workers_extract: Maximum parallel workers for extraction (network I/O intensive).
        :param max_workers_shuffle: Maximum parallel workers for shuffling (CPU/memory intensive).
        :param cleanup_temp: Whether to delete temp directory after shuffling.
        :param verbose: If False, suppress INFO level logging in parallel workers.

        :raise SomaExtractError: If SOMA reads fail
        :raise IOError: If file operations fail
        :raise ValueError: If output_format is invalid
        """
        _extract.run_extract_randomized(
            curriculum_metadata=curriculum_metadata,
            output_dir=output_dir,
            partition_index=partition_index,
            max_ranges_per_partition=max_ranges_per_partition,
            output_format=output_format,
            temp_dir=temp_dir,
            max_workers_extract=max_workers_extract,
            max_workers_shuffle=max_workers_shuffle,
            cleanup_temp=cleanup_temp,
            verbose=verbose,
        )

    def extract_grouped(
        self,
        *,
        curriculum_metadata: GroupedCurriculumMetadata,
        output_dir: Path,
        partition_index: int = 0,
        max_bins_per_partition: int | None = None,
        output_format: Literal["zarr", "h5ad"] = "h5ad",
        max_workers: int | None = None,
        verbose: bool = False,
    ) -> None:
        """
        Extract grouped bins directly to output files.

        Support distributed execution via partition_index and max_bins_per_partition.
        No shuffle stage — cells from the same group stay together.

        :param curriculum_metadata: Metadata with grouped_bins
        :param output_dir: Output directory for extract files
        :param partition_index: Zero-based partition index for distributed workers
        :param max_bins_per_partition: Number of bins per partition (for distribution)
        :param output_format: Output format - "zarr" or "h5ad"
        :param max_workers: Parallel workers within this partition
        :param verbose: Enable verbose logging in workers

        :raise SomaExtractError: If SOMA reads fail
        :raise IOError: If file operations fail
        :raise ValueError: If output_format is invalid or grouped_bins is None
        """
        _extract.extract_grouped_bins(
            curriculum_metadata=curriculum_metadata,
            output_dir=output_dir,
            partition_index=partition_index,
            max_bins_per_partition=max_bins_per_partition,
            output_format=output_format,
            max_workers=max_workers,
            verbose=verbose,
        )
