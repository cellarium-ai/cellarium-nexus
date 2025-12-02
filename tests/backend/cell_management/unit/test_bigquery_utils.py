import types

import pytest

from cellarium.nexus.backend.cell_management.utils import bigquery_utils


def test_canonical_json_is_stable_and_sorted() -> None:
    """
    Assert that _canonical_json produces deterministic, compact JSON with sorted keys.
    """
    a = {"b": 1, "a": 2, "nested": {"y": 1, "x": 2}}
    b = {"nested": {"x": 2, "y": 1}, "a": 2, "b": 1}

    ja = bigquery_utils._canonical_json(data=a)
    jb = bigquery_utils._canonical_json(data=b)

    assert ja == jb
    # Should be compact (no spaces after commas/colons) and sorted by keys
    assert ":" in ja and "," in ja
    assert " " not in ja


def test_hash_to_key_prefix_and_md5_length() -> None:
    """
    Assert that _hash_to_key prefixes with namespace and uses 32-char md5 digest.
    """
    key = bigquery_utils._hash_to_key(prefix="countcache", seed="seed")
    assert key.startswith("countcache:")
    digest = key.split(":", 1)[1]
    assert len(digest) == 32
    # hex chars
    int(digest, 16)


@pytest.fixture()
def cache_namespace() -> str:
    """
    Provide a test cache namespace for manager construction.

    :return: Cache namespace string
    """
    return "test-namespace"


class MockDataOperator:
    """
    Mock data operator implementing DataOperatorProtocol for testing.
    """

    def __init__(self, *, impls: dict | None = None, recorder: types.SimpleNamespace | None = None) -> None:
        self.impls = impls or {}
        self.recorder = recorder or types.SimpleNamespace(
            count_cells=0, get_categorical_obs_columns=0, get_distinct_obs_values=0
        )

    def count_cells(self, *, filter_statements: dict | None = None) -> int:
        """Return a value and record call count."""
        self.recorder.count_cells += 1
        fn = self.impls.get("count_cells")
        return 123 if fn is None else fn(filter_statements=filter_statements)

    def get_categorical_obs_columns(self, *, threshold: int, exclude: list[str] | None = None) -> set[str]:
        """Return a set and record call count."""
        self.recorder.get_categorical_obs_columns += 1
        fn = self.impls.get("get_categorical_obs_columns")
        return {"organism"} if fn is None else fn(threshold=threshold, exclude=exclude)

    def get_distinct_obs_values(self, *, column_name: str) -> list[str]:
        """Return a list and record call count."""
        self.recorder.get_distinct_obs_values += 1
        fn = self.impls.get("get_distinct_obs_values")
        return ["Homo sapiens"] if fn is None else fn(column_name=column_name)


def test_get_cached_count_uses_cache(cache_namespace: str) -> None:
    """
    Assert first call hits operator and subsequent calls read cache for same inputs.
    Also assert canonical JSON means dict order does not affect the cache key.
    """
    rec = types.SimpleNamespace(count_cells=0, get_categorical_obs_columns=0, get_distinct_obs_values=0)
    operator = MockDataOperator(
        impls={"count_cells": lambda **kwargs: 42},
        recorder=rec,
    )

    mgr = bigquery_utils.OmicsCachedDataManager(operator=operator, cache_namespace=cache_namespace)

    filters_a = {"b": 1, "a": 2}
    filters_b = {"a": 2, "b": 1}  # same content, different order

    # First call should invoke operator
    out1 = mgr.get_cached_count(filters_dict=filters_a)
    assert out1 == 42
    assert rec.count_cells == 1

    # Second call with same inputs should use cache (no additional operator call)
    out2 = mgr.get_cached_count(filters_dict=filters_a)
    assert out2 == 42
    assert rec.count_cells == 1

    # Same filter content, different order should still be cached
    out3 = mgr.get_cached_count(filters_dict=filters_b)
    assert out3 == 42
    assert rec.count_cells == 1


def test_get_cached_categorical_obs_columns_caches(cache_namespace: str) -> None:
    """
    Assert categorical columns are cached per namespace/threshold.
    """
    rec = types.SimpleNamespace(count_cells=0, get_categorical_obs_columns=0, get_distinct_obs_values=0)
    operator = MockDataOperator(
        impls={"get_categorical_obs_columns": lambda **kwargs: {"organism", "tissue"}},
        recorder=rec,
    )

    mgr = bigquery_utils.OmicsCachedDataManager(operator=operator, cache_namespace=cache_namespace)

    out1 = mgr.get_cached_categorical_obs_columns(distinct_threshold=100)
    assert out1 == {"organism", "tissue"}
    assert rec.get_categorical_obs_columns == 1

    # Repeat same request → served from cache
    out2 = mgr.get_cached_categorical_obs_columns(distinct_threshold=100)
    assert out2 == {"organism", "tissue"}
    assert rec.get_categorical_obs_columns == 1

    # Change threshold → cache miss
    out3 = mgr.get_cached_categorical_obs_columns(distinct_threshold=200)
    assert out3 == {"organism", "tissue"}
    assert rec.get_categorical_obs_columns == 2


def test_get_cached_distinct_obs_values_caches(cache_namespace: str) -> None:
    """
    Assert distinct values are cached per namespace/column.
    """
    rec = types.SimpleNamespace(count_cells=0, get_categorical_obs_columns=0, get_distinct_obs_values=0)
    operator = MockDataOperator(
        impls={"get_distinct_obs_values": lambda **kwargs: ["A", "B"]},
        recorder=rec,
    )

    mgr = bigquery_utils.OmicsCachedDataManager(operator=operator, cache_namespace=cache_namespace)

    out1 = mgr.get_cached_distinct_obs_values(column_name="organism")
    assert out1 == ["A", "B"]
    assert rec.get_distinct_obs_values == 1

    # Repeat same request → served from cache
    out2 = mgr.get_cached_distinct_obs_values(column_name="organism")
    assert out2 == ["A", "B"]
    assert rec.get_distinct_obs_values == 1

    # Different column → cache miss
    out3 = mgr.get_cached_distinct_obs_values(column_name="tissue")
    assert out3 == ["A", "B"]
    assert rec.get_distinct_obs_values == 2
