"""
SOMA grouped curriculum metadata utilities.

This module provides planning logic for grouped SOMA extracts where cells
from the same group stay together.
"""

import logging

import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops import exceptions
from cellarium.nexus.omics_datastore.soma_ops.curriculum_randomized import read_var_joinids
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter
from cellarium.nexus.shared.schemas.omics_datastore import GroupedBin, SomaCurriculumMetadata

logger = logging.getLogger(__name__)


def _build_group_filter(group_columns: list[str], group_values: tuple) -> str:
    """
    Build a SOMA value_filter expression for a specific group.

    :param group_columns: List of column names used for grouping
    :param group_values: Tuple of values corresponding to each column

    :return: SOMA value_filter expression string
    """
    conditions = []
    for col, val in zip(group_columns, group_values):
        # Escape quotes in string values
        if isinstance(val, str):
            escaped_val = val.replace('"', '\\"')
            conditions.append(f'{col} == "{escaped_val}"')
        else:
            conditions.append(f"{col} == {val}")

    return " and ".join(conditions)


def prepare_grouped_curriculum(
    *,
    experiment_uri: str,
    filters: dict[str, object] | None,
    extract_bin_keys: list[str],
    bin_size: int,
    var_filter_column: str | None = None,
    var_filter_values: list[str] | None = None,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    x_layer: str = "X",
) -> SomaCurriculumMetadata:
    """
    Prepare curriculum metadata for grouped extraction.

    Compute grouped bins where cells from the same group stay together,
    with each bin not exceeding bin_size.

    :param experiment_uri: URI of the SOMA experiment
    :param filters: Dict of filter conditions using the Nexus format
    :param extract_bin_keys: List of obs column names to group by
    :param bin_size: Maximum number of cells per bin
    :param var_filter_column: Name of the var column to filter features by
    :param var_filter_values: List of values to match in the var filter column
    :param obs_columns: List of obs columns to include in extraction
    :param var_columns: List of var columns to include in extraction
    :param x_layer: Name of the SOMA X layer to read counts from

    :raise SomaReadError: If SOMA reads fail
    :raise SomaFilterError: If filter translation fails
    :raise ValueError: If bin_size is not positive or extract_bin_keys is empty

    :return: SomaCurriculumMetadata with grouped_bins
    """
    if bin_size <= 0:
        raise ValueError(f"bin_size must be positive, got {bin_size}")
    if not extract_bin_keys:
        raise ValueError("extract_bin_keys cannot be empty")

    value_filter = build_soma_value_filter(filters=filters)

    logger.info(
        f"Preparing grouped curriculum for {experiment_uri} with bin_size={bin_size},\n"
        f"extract_bin_keys={extract_bin_keys},\n"
        f"filter={value_filter if value_filter else 'no filter'}"
    )

    # Compute var_ids if feature filtering is requested
    var_ids: list[int] | None = None
    if var_filter_column and var_filter_values:
        var_ids = read_var_joinids(
            experiment_uri=experiment_uri,
            var_filter_column=var_filter_column,
            var_filter_values=var_filter_values,
        )
        logger.info(f"Found {len(var_ids)} var joinids for feature filter")

    # Read obs with grouping columns
    try:
        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            columns_to_read = ["soma_joinid"] + extract_bin_keys
            obs_query = exp.obs.read(
                column_names=columns_to_read,
                value_filter=value_filter if value_filter else None,
            )
            obs_df = obs_query.concat().to_pandas()
    except Exception as e:
        raise exceptions.SomaReadError(f"Failed to read obs from {experiment_uri}: {e}") from e

    total_cells = len(obs_df)
    logger.info(f"Found {total_cells} cells to extract")

    if total_cells == 0:
        raise exceptions.SomaPrepareCurriculumMetadataError(f"No cells found matching the filter: {value_filter}")

    # Group by extract_bin_keys and create bins
    grouped_bins: list[GroupedBin] = []

    for group_values, group_df in obs_df.groupby(extract_bin_keys, sort=True):
        # Handle single column grouping (returns scalar, not tuple)
        if not isinstance(group_values, tuple):
            group_values = (group_values,)

        group_filter = _build_group_filter(extract_bin_keys, group_values)
        group_key = "||".join(str(v) for v in group_values)

        # Sort joinids within group for deterministic chunking
        joinids = group_df["soma_joinid"].sort_values().values

        # Split large groups into multiple bins
        for chunk_start in range(0, len(joinids), bin_size):
            chunk_joinids = joinids[chunk_start : chunk_start + bin_size]
            grouped_bins.append(
                GroupedBin(
                    group_key=group_key,
                    group_filter=group_filter,
                    joinid_min=int(chunk_joinids.min()),
                    joinid_max=int(chunk_joinids.max()),
                    cell_count=len(chunk_joinids),
                )
            )

    num_grouped_bins = len(grouped_bins)
    logger.info(f"Created {num_grouped_bins} grouped bins")

    metadata = SomaCurriculumMetadata(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        total_cells=total_cells,
        filters=filters,
        var_joinids=var_ids,
        var_filter_column=var_filter_column,
        var_filter_values=var_filter_values,
        obs_columns=obs_columns,
        var_columns=var_columns,
        x_layer=x_layer,
        extract_bin_keys=extract_bin_keys,
        grouped_bins=grouped_bins,
        num_grouped_bins=num_grouped_bins,
    )

    logger.info(
        f"Grouped curriculum metadata created: {metadata.total_cells} cells " f"in {metadata.num_grouped_bins} bins"
    )

    return metadata
