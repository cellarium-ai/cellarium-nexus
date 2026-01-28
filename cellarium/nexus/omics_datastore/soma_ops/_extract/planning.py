"""
SOMA curriculum planning utilities.

This module provides high-level planning logic for SOMA extracts, dispatching to
specific planning strategies (randomized or grouped).
"""

from typing import Any

from cellarium.nexus.omics_datastore.soma_ops._extract.prepare_curriculum_grouped import prepare_grouped_curriculum
from cellarium.nexus.omics_datastore.soma_ops._extract.prepare_curriculum_randomized import prepare_extract_curriculum
from cellarium.nexus.shared.schemas.omics_datastore import (
    GroupedCurriculumMetadata,
    RandomizedCurriculumMetadata,
)


def plan_curriculum(
    *,
    experiment_uri: str,
    filters: dict[str, Any] | None,
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

    Dispatches to grouped or randomized curriculum planning based on extract_bin_keys.

    :param experiment_uri: URI of the SOMA experiment.
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
    if extract_bin_keys:
        # Grouped extraction mode
        if bin_size is None:
            raise ValueError("bin_size is required for grouped extraction")
        return prepare_grouped_curriculum(
            experiment_uri=experiment_uri,
            filters=filters,
            extract_bin_keys=extract_bin_keys,
            bin_size=bin_size,
            var_filter_column=var_filter_column,
            var_filter_values=var_filter_values,
            obs_columns=obs_columns,
            var_columns=var_columns,
            x_layer=x_layer,
        )
    else:
        # Randomized extraction mode
        if range_size is None:
            raise ValueError("range_size is required for randomized extraction")
        if extract_bin_size is None:
            raise ValueError("extract_bin_size is required for randomized extraction")
        return prepare_extract_curriculum(
            experiment_uri=experiment_uri,
            filters=filters,
            range_size=range_size,
            extract_bin_size=extract_bin_size,
            shuffle_ranges=shuffle_ranges,
            var_filter_column=var_filter_column,
            var_filter_values=var_filter_values,
            obs_columns=obs_columns,
            var_columns=var_columns,
            x_layer=x_layer,
        )
