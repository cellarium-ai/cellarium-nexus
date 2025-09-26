from __future__ import annotations

import pytest

from cellarium.nexus.backend.cell_management.utils import bigquery_utils


class BigQueryCachedManagerStub:
    """
    Provide deterministic responses for `BigQueryCachedDataManager` methods.

    :param args: Positional constructor arguments (ignored)
    :param kwargs: Keyword constructor arguments (ignored)
    """

    def __init__(self, *args, **kwargs) -> None:  # noqa: D401
        self.args = args
        self.kwargs = kwargs

    def get_cached_count_bq(
        self,
        *,
        dataset_name: str,
        filters_dict: dict[str, object],
        timeout: int = 0,
    ) -> int:
        return 42

    def get_cached_categorical_columns_bq(
        self,
        *,
        dataset_name: str,
        table_name: str = "",
        distinct_threshold: int = 0,
        timeout: int = 0,
    ) -> set[str]:
        return set()

    def get_cached_distinct_values_bq(
        self,
        *,
        dataset_name: str,
        column_name: str,
        table_name: str = "",
        limit: int = 0,
        timeout: int = 0,
    ) -> list[str]:
        return []


@pytest.fixture()
def bigquery_cached_manager_stub(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Patch `BigQueryCachedDataManager` to a deterministic stub for backend tests.

    :param monkeypatch: Pytest monkeypatch fixture
    """

    monkeypatch.setattr(bigquery_utils, "BigQueryCachedDataManager", BigQueryCachedManagerStub)
