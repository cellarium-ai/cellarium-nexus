"""
SOMA query utilities.

This module provides low-level query operations for SOMA experiments.
"""

import logging
from typing import Any

import pyarrow as pa
import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaReadError
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter

logger = logging.getLogger(__name__)


def count_obs(*, experiment_uri: str, filter_statements: dict[str, Any] | None = None) -> int:
    """
    Count cells in the SOMA obs table matching the given filters.

    :param experiment_uri: URI of the SOMA experiment
    :param filter_statements: Dict of filter conditions using the Nexus format

    :raise SomaReadError: If SOMA read fails

    :return: Total number of matching cells
    """
    try:
        value_filter = build_soma_value_filter(filters=filter_statements)

        logger.info(f"Counting cells with filter: {value_filter if value_filter else 'no filter'}")

        with tiledbsoma.open(experiment_uri, mode="r") as exp:
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


def get_obs_string_columns(*, experiment_uri: str, exclude: list[str] | None = None) -> set[str]:
    """
    Determine categorical string columns in obs data.

    Identifies plain strings, large strings, or dictionary-encoded strings.

    :param experiment_uri: URI of the SOMA experiment
    :param exclude: Optional list of column names to skip

    :raise SomaReadError: If SOMA read fails

    :return: Set of column names that are string-like
    """
    try:
        exclude_set = set(exclude or [])

        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            obs_schema = exp.obs.schema

            # Get string-like columns from schema (plain, large, or dictionary-encoded strings)
            string_columns = set()
            for field in obs_schema:
                if field.name in exclude_set:
                    continue

                ftype = field.type
                is_string_like = bool(
                    pa.types.is_string(ftype)
                    # or pa.types.is_large_string(ftype)
                    or (pa.types.is_dictionary(ftype) and pa.types.is_string(ftype.value_type))
                )

                logger.debug(f"Field {field.name}: type={ftype}, is_string_like={is_string_like}")

                if is_string_like:
                    string_columns.add(field.name)

            logger.info(f"Found {len(string_columns)} string-like columns in SOMA obs: {string_columns}")

            return string_columns

    except Exception as e:
        logger.error(f"Failed to get categorical obs columns: {e}")
        raise SomaReadError("SOMA categorical columns operation failed") from e


def get_obs_distinct_values(*, experiment_uri: str, column_name: str) -> list[str]:
    """
    Fetch distinct values for an obs column.

    :param experiment_uri: URI of the SOMA experiment
    :param column_name: Column to retrieve values for

    :raise SomaReadError: If SOMA read fails

    :return: List of distinct values
    """
    try:
        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            obs_query = exp.obs.read(column_names=[column_name])
            obs_df = obs_query.concat().to_pandas()

            # Get unique values, drop nulls, convert to string
            unique_values = obs_df[column_name].dropna().unique()
            result = [str(v) for v in unique_values]

            return result

    except Exception as e:
        logger.error(f"Failed to get distinct obs values: {e}")
        raise SomaReadError("SOMA distinct values operation failed") from e
