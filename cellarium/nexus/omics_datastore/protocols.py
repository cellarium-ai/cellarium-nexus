"""
Protocol definitions for omics datastore operators.
"""

from typing import Any, Protocol


class DataOperatorProtocol(Protocol):
    """
    Define required methods for data operators used by OmicsCachedDataManager.

    Both BigQueryDataOperator and TileDBSOMADataOperator implement this protocol
    to enable interchangeable use with the cached data manager.
    """

    def count_cells(self, *, filter_statements: dict[str, Any] | None = None) -> int:
        """
        Count cells matching the given filters.

        :param filter_statements: Optional dictionary of filters to apply

        :raise Exception: If query execution fails

        :return: Total count of matching cells
        """
        ...

    def get_categorical_obs_columns(
        self,
        *,
        threshold: int,
        exclude: list[str] | None = None,
    ) -> set[str]:
        """
        Determine categorical string columns in obs/cell_info by distinct count threshold.

        :param threshold: Maximum distinct values to consider a column categorical
        :param exclude: Optional list of column names to skip

        :raise Exception: If query execution fails

        :return: Set of column names that are categorical
        """
        ...

    def get_distinct_obs_values(
        self,
        *,
        column_name: str,
    ) -> list[str]:
        """
        Fetch distinct values for an obs/cell_info column.

        :param column_name: Column to retrieve values for

        :raise Exception: If query execution fails

        :return: List of distinct values
        """
        ...
