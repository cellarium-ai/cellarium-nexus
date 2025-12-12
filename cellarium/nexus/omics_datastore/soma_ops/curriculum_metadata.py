"""
SOMA extract planning utilities.

This module provides high-level planning logic for SOMA extracts.
"""

import logging
import math
import random
import typing as t

import numpy as np
import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops import exceptions
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter
from cellarium.nexus.shared.schemas.omics_datastore import IdContiguousRange, SomaCurriculumMetadata

logger = logging.getLogger(__name__)

T = t.TypeVar("T")


def read_filtered_joinids(
    *,
    experiment_uri: str,
    value_filter: str,
) -> np.ndarray:
    """
    Read and sort soma_joinids matching the given filter.

    :param experiment_uri: URI of the SOMA experiment
    :param value_filter: SOMA obs value_filter expression

    :raise SomaReadError: If SOMA read fails

    :return: Sorted array of soma_joinids
    """
    try:
        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            # Read soma_joinid column with the filter
            obs_query = exp.obs.read(
                column_names=["soma_joinid"],
                value_filter=value_filter if value_filter else None,
            )

            # Get joinids as numpy array
            obs_df = obs_query.concat().to_pandas()
            joinids = obs_df["soma_joinid"].to_numpy(dtype=np.int64)

            # Sort joinids
            return np.sort(joinids)
    except Exception as e:
        raise exceptions.SomaReadError(f"Failed to read joinids from {experiment_uri}: {e}") from e


def shuffle_but_keep_last(*, values: list[T]) -> list[T]:
    """
    Shuffle all values in the list except for the last one, the last value remains the same as it used to be

    :param values: Input list of values to be shuffled.

    :return: List of shuffled values, where the last member is the same as it was.
    """
    # keep last element fixed
    last = values[-1]

    # shuffle everything except last
    to_shuffle = values[:-1]
    shuffled = random.sample(to_shuffle, k=len(to_shuffle))
    shuffled.append(last)
    return shuffled


def compute_contiguous_ranges(*, values: np.ndarray, chunk_size: int, shuffle: bool = False) -> list[IdContiguousRange]:
    """
    Compute contiguous ranges from sorted values.

    Split the sorted values into chunks of the specified size and create
    a range for each chunk, where each range captures the first and last
    value in that chunk.

    :param values: Sorted array of integer values
    :param chunk_size: Target number of values per range
    :param shuffle: Whether to shuffle ranges (all except last)

    :raise ValueError: If chunk_size is not positive

    :return: Tuple of (list of ranges, index of incomplete range or None)
    """
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    total_count = len(values)
    if total_count == 0:
        return []

    ranges = []

    for i in range(0, total_count, chunk_size):
        chunk = values[i : i + chunk_size]
        value_range = IdContiguousRange(
            start=int(chunk[0]),
            end=int(chunk[-1]),
        )
        ranges.append(value_range)

    if shuffle:
        ranges = shuffle_but_keep_last(values=ranges)

    return ranges


def compute_output_ids(*, n_values: int, chunk_size: int, shuffle: bool = False) -> list[int]:
    num_output_chunks = math.ceil(n_values / chunk_size)
    output_ids = list(range(num_output_chunks))

    if shuffle:
        output_ids = shuffle_but_keep_last(values=output_ids)

    return output_ids


def read_var_joinids(
    *,
    experiment_uri: str,
    var_filter_column: str,
    var_filter_values: list[str],
) -> list[int]:
    """
    Read var soma_joinids matching the given filter values.

    Return joinids in the same order as var_filter_values to preserve
    the requested feature ordering.

    :param experiment_uri: URI of the SOMA experiment
    :param var_filter_column: Name of the var column to filter on
    :param var_filter_values: List of values to match in the filter column

    :raise SomaReadError: If SOMA read fails

    :return: List of soma_joinids in the order of var_filter_values
    """
    try:
        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            rna = exp.ms["RNA"]
            var_df = rna.var.read(column_names=["soma_joinid", var_filter_column]).concat().to_pandas()

            # Create mapping from filter value to joinid
            value_to_joinid = dict(zip(var_df[var_filter_column], var_df["soma_joinid"]))

            # Return joinids in the order of var_filter_values
            result = []
            for value in var_filter_values:
                if value in value_to_joinid:
                    result.append(int(value_to_joinid[value]))
                else:
                    logger.warning(f"Feature value '{value}' not found in var column '{var_filter_column}'")

            return result
    except Exception as e:
        raise exceptions.SomaReadError(f"Failed to read var joinids from {experiment_uri}: {e}") from e


def prepare_extract_curriculum(
    *,
    experiment_uri: str,
    filters: dict[str, object] | None,
    range_size: int,
    output_chunk_size: int,
    shuffle_ranges: bool = True,
    shuffle_output_ids: bool = True,
    var_filter_column: str | None = None,
    var_filter_values: list[str] | None = None,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    x_layer: str = "X",
) -> SomaCurriculumMetadata:
    """
    Plan a SOMA extract by computing soma_joinid ranges for a given filter.

    Read filtered obs data, compute contiguous joinid ranges, and return
    an extract metadata. This is the main entry point for dashboard/backend
    to metadata SOMA extracts.

    :param experiment_uri: URI of the SOMA experiment.
    :param filters: Dict of filter conditions using the Nexus format.
    :param range_size: Target number of cells per range (for extraction).
    :param output_chunk_size: Target number of cells per output chunk (for shuffling).
    :param shuffle_ranges: Whether to shuffle id ranges.
    :param shuffle_output_ids: Whether to shuffle output chunk ids.
    :param var_filter_column: Name of the var column to filter features by.
    :param var_filter_values: List of values to match in the var filter column.
    :param obs_columns: List of obs columns to include in extraction.
    :param var_columns: List of var columns to include in extraction.
    :param x_layer: Name of the SOMA X layer to read counts from.

    :raise SomaReadError: If SOMA reads fail
    :raise SomaFilterError: If filter translation fails
    :raise ValueError: If range_size or output_chunk_size is not positive

    :return: SomaCurriculumMetadata object
    """
    if range_size <= 0:
        raise ValueError(f"range_size must be positive, got {range_size}")
    if output_chunk_size <= 0:
        raise ValueError(f"output_chunk_size must be positive, got {output_chunk_size}")

    value_filter = build_soma_value_filter(filters=filters)

    logger.info(
        f"Planning SOMA extract for {experiment_uri} with range_size={range_size},\n"
        f"filter={value_filter if value_filter else 'no filter'},\n"
        f"var_filter_column={var_filter_column},\n"
        f"var_filter_values_count={len(var_filter_values) if var_filter_values else 0}"
    )

    # Compute var_ids if feature filtering is requested
    var_ids: list[int] | None = None
    if var_filter_column and var_filter_values:
        var_ids = read_var_joinids(
            experiment_uri=experiment_uri,
            var_filter_column=var_filter_column,
            var_filter_values=var_filter_values,
        )
        logger.info(f"Found {len(var_ids)} var obs_ids for feature filter")

    # Read and sort filtered obs_ids
    obs_ids = read_filtered_joinids(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
    )

    total_cells = len(obs_ids)
    logger.info(f"Found {total_cells} cells to extract")

    if total_cells == 0:
        raise exceptions.SomaPrepareCurriculumMetadataError(f"No cells found matching the filter: {value_filter}")

    # Compute contiguous id ranges
    id_ranges = compute_contiguous_ranges(values=obs_ids, chunk_size=range_size, shuffle=shuffle_ranges)
    output_ids = compute_output_ids(n_values=total_cells, chunk_size=output_chunk_size, shuffle=shuffle_output_ids)
    logger.info(f"Computed {len(id_ranges)} ranges. Output chunk n: {output_ids}")

    num_output_chunks = len(output_ids)
    num_ranges = len(id_ranges)
    last_chunk_size = total_cells % output_chunk_size

    metadata = SomaCurriculumMetadata(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        id_ranges=id_ranges,
        total_cells=total_cells,
        range_size=range_size,
        num_ranges=num_ranges,
        output_chunk_size=output_chunk_size,
        num_output_chunks=num_output_chunks,
        last_chunk_size=last_chunk_size,
        output_chunk_indexes=output_ids,
        filters=filters,
        var_joinids=var_ids,
        var_filter_column=var_filter_column,
        var_filter_values=var_filter_values,
        obs_columns=obs_columns,
        var_columns=var_columns,
        x_layer=x_layer,
    )

    logger.info(f"Extract metadata created: {metadata.total_cells} cells in {len(metadata.id_ranges)} ranges")

    return metadata
