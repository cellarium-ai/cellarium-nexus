"""
SOMA extract planning utilities.

This module provides high-level planning logic for SOMA extracts.
"""

import logging
import random

import numpy as np
import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaReadError
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter
from cellarium.nexus.shared.schemas.omics_datastore import SomaExtractPlan, SomaJoinIdRange

logger = logging.getLogger(__name__)


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
        raise SomaReadError(f"Failed to read joinids from {experiment_uri}: {e}") from e


def compute_contiguous_ranges(
    *,
    values: np.ndarray,
    chunk_size: int,
) -> list[SomaJoinIdRange]:
    """
    Compute contiguous ranges from sorted values.

    Split the sorted values into chunks of the specified size and create
    a range for each chunk, where each range captures the first and last
    value in that chunk.

    :param values: Sorted array of integer values
    :param chunk_size: Target number of values per range

    :raise ValueError: If chunk_size is not positive

    :return: List of ranges, each containing start and end values
    """
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    total_count = len(values)
    if total_count == 0:
        return []

    ranges = []
    for i in range(0, total_count, chunk_size):
        chunk = values[i : i + chunk_size]
        value_range = SomaJoinIdRange(
            start=int(chunk[0]),
            end=int(chunk[-1]),
        )
        ranges.append(value_range)

    return ranges


def plan_soma_extract(
    *,
    experiment_uri: str,
    filters: dict[str, object] | None,
    range_size: int,
    output_chunk_size: int,
    shuffle_ranges: bool = True,
) -> SomaExtractPlan:
    """
    Plan a SOMA extract by computing soma_joinid ranges for a given filter.

    Read filtered obs data, compute contiguous joinid ranges, and return
    an extract plan. This is the main entry point for dashboard/backend
    to plan SOMA extracts.

    :param experiment_uri: URI of the SOMA experiment
    :param filters: Dict of filter conditions using the Nexus format
    :param range_size: Target number of cells per range (for extraction)
    :param output_chunk_size: Target number of cells per output chunk (for shuffling)
    :param shuffle_ranges: Whether to shuffle joinid ranges

    :raise SomaReadError: If SOMA reads fail
    :raise SomaFilterError: If filter translation fails
    :raise ValueError: If range_size or output_chunk_size is not positive

    :return: SomaExtractPlan object
    """
    if range_size <= 0:
        raise ValueError(f"range_size must be positive, got {range_size}")
    if output_chunk_size <= 0:
        raise ValueError(f"output_chunk_size must be positive, got {output_chunk_size}")

    value_filter = build_soma_value_filter(filters=filters)

    logger.info(
        f"Planning SOMA extract for {experiment_uri} with range_size={range_size}, "
        f"filter={value_filter if value_filter else 'no filter'}"
    )

    # Read and sort filtered joinids
    joinids = read_filtered_joinids(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
    )

    total_cells = len(joinids)
    logger.info(f"Found {total_cells} cells to extract")

    if total_cells == 0:
        logger.warning("No cells found matching the filter")
        return SomaExtractPlan(
            experiment_uri=experiment_uri,
            value_filter=value_filter,
            joinid_ranges=[],
            total_cells=0,
            range_size=range_size,
            output_chunk_size=output_chunk_size,
            filters=filters,
        )

    # Compute joinid ranges
    joinid_ranges = compute_contiguous_ranges(
        values=joinids,
        chunk_size=range_size,
    )

    logger.info(f"Created {len(joinid_ranges)} joinid ranges")

    # Shuffle ranges if requested
    if shuffle_ranges:
        random.shuffle(joinid_ranges)
        logger.info("Shuffled joinid ranges")

    plan = SomaExtractPlan(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        joinid_ranges=joinid_ranges,
        total_cells=total_cells,
        range_size=range_size,
        output_chunk_size=output_chunk_size,
        filters=filters,
    )

    logger.info(f"Extract plan created: {plan.total_cells} cells in {len(plan.joinid_ranges)} ranges")

    return plan
