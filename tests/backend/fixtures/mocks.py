from __future__ import annotations

import pytest

from cellarium.nexus.backend.cell_management.utils import bigquery_utils


class OmicsCachedManagerStub:
    """
    Provide deterministic responses for `OmicsCachedDataManager` methods.

    :param operator: Data operator (ignored in stub)
    :param cache_namespace: Cache namespace (ignored in stub)
    """

    def __init__(self, *, operator: object = None, cache_namespace: str = "") -> None:
        self.operator = operator
        self.cache_namespace = cache_namespace

    def get_cached_count(
        self,
        *,
        filters_dict: dict[str, object],
        timeout: int = 0,
    ) -> int:
        return 42

    def get_cached_categorical_obs_columns(
        self,
        *,
        distinct_threshold: int = 0,
        timeout: int = 0,
    ) -> set[str]:
        return set()

    def get_cached_distinct_obs_values(
        self,
        *,
        column_name: str,
        timeout: int = 0,
    ) -> list[str]:
        return []


@pytest.fixture()
def omics_cached_manager_stub(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Patch `OmicsCachedDataManager` to a deterministic stub for backend tests.

    :param monkeypatch: Pytest monkeypatch fixture
    """

    monkeypatch.setattr(bigquery_utils, "OmicsCachedDataManager", OmicsCachedManagerStub)
