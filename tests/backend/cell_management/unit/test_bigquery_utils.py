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
def project_id() -> str:
    """
    Provide a test project id for manager construction.

    :return: Project ID string
    """
    return "test-project"


def _patch_bq_operator(monkeypatch: pytest.MonkeyPatch, **impls):
    """
    Patch BigQueryDataOperator methods used by BigQueryCachedDataManager.

    :param monkeypatch: Pytest monkeypatch fixture

    :return: A recorder namespace with call counters
    """
    # Recorder object to track calls
    rec = types.SimpleNamespace(count_cells=0, get_categorical_string_columns=0, get_distinct_values=0)

    class _Op:
        def __init__(self, *args, **kwargs) -> None:  # noqa: D401
            """Accept any args without side effects."""

        def count_cells(self, *args, **kwargs) -> int:  # noqa: D401
            """Return a value and record call count."""
            rec.count_cells += 1
            fn = impls.get("count_cells")
            return 123 if fn is None else fn(**kwargs)

        def get_categorical_string_columns(self, *args, **kwargs) -> set[str]:  # noqa: D401
            """Return a set and record call count."""
            rec.get_categorical_string_columns += 1
            fn = impls.get("get_categorical_string_columns")
            return {"organism"} if fn is None else fn(**kwargs)

        def get_distinct_values(self, *args, **kwargs) -> list[str]:  # noqa: D401
            """Return a list and record call count."""
            rec.get_distinct_values += 1
            fn = impls.get("get_distinct_values")
            return ["Homo sapiens"] if fn is None else fn(**kwargs)

    monkeypatch.setattr(bigquery_utils.bq_ops, "BigQueryDataOperator", _Op, raising=True)
    return rec


def test_get_cached_count_bq_uses_cache(monkeypatch: pytest.MonkeyPatch, project_id: str) -> None:
    """
    Assert first call hits operator and subsequent calls read cache for same inputs.
    Also assert canonical JSON means dict order does not affect the cache key.
    """
    rec = _patch_bq_operator(monkeypatch, count_cells=lambda **kwargs: 42)

    mgr = bigquery_utils.BigQueryCachedDataManager(project_id=project_id)

    filters_a = {"b": 1, "a": 2}
    filters_b = {"a": 2, "b": 1}  # same content, different order

    # First call should invoke operator
    out1 = mgr.get_cached_count_bq(dataset_name="ds1", filters_dict=filters_a)
    assert out1 == 42
    assert rec.count_cells == 1

    # Second call with same inputs should use cache (no additional operator call)
    out2 = mgr.get_cached_count_bq(dataset_name="ds1", filters_dict=filters_a)
    assert out2 == 42
    assert rec.count_cells == 1

    # Same filter content, different order should still be cached
    out3 = mgr.get_cached_count_bq(dataset_name="ds1", filters_dict=filters_b)
    assert out3 == 42
    assert rec.count_cells == 1

    # Different dataset should miss cache and call operator
    out4 = mgr.get_cached_count_bq(dataset_name="ds2", filters_dict=filters_a)
    assert out4 == 42
    assert rec.count_cells == 2


def test_get_cached_categorical_columns_bq_caches(monkeypatch: pytest.MonkeyPatch, project_id: str) -> None:
    """
    Assert categorical columns are cached per dataset/table/threshold namespace.
    """
    rec = _patch_bq_operator(
        monkeypatch,
        get_categorical_string_columns=lambda **kwargs: {"organism", "tissue"},
    )

    mgr = bigquery_utils.BigQueryCachedDataManager(project_id=project_id)

    out1 = mgr.get_cached_categorical_columns_bq(dataset_name="ds1", table_name="cell_info", distinct_threshold=100)
    assert out1 == {"organism", "tissue"}
    assert rec.get_categorical_string_columns == 1

    # Repeat same request → served from cache
    out2 = mgr.get_cached_categorical_columns_bq(dataset_name="ds1", table_name="cell_info", distinct_threshold=100)
    assert out2 == {"organism", "tissue"}
    assert rec.get_categorical_string_columns == 1

    # Change threshold → cache miss
    out3 = mgr.get_cached_categorical_columns_bq(dataset_name="ds1", table_name="cell_info", distinct_threshold=200)
    assert out3 == {"organism", "tissue"}
    assert rec.get_categorical_string_columns == 2


def test_get_cached_distinct_values_bq_caches(monkeypatch: pytest.MonkeyPatch, project_id: str) -> None:
    """
    Assert distinct values are cached per dataset/table/column/limit.
    """
    rec = _patch_bq_operator(
        monkeypatch,
        get_distinct_values=lambda **kwargs: ["A", "B"],
    )

    mgr = bigquery_utils.BigQueryCachedDataManager(project_id=project_id)

    out1 = mgr.get_cached_distinct_values_bq(
        dataset_name="ds1",
        column_name="organism",
        table_name="cell_info",
        limit=100,
    )
    assert out1 == ["A", "B"]
    assert rec.get_distinct_values == 1

    # Repeat same request → served from cache
    out2 = mgr.get_cached_distinct_values_bq(
        dataset_name="ds1",
        column_name="organism",
        table_name="cell_info",
        limit=100,
    )
    assert out2 == ["A", "B"]
    assert rec.get_distinct_values == 1

    # Change limit → cache miss
    out3 = mgr.get_cached_distinct_values_bq(
        dataset_name="ds1",
        column_name="organism",
        table_name="cell_info",
        limit=200,
    )
    assert out3 == ["A", "B"]
    assert rec.get_distinct_values == 2
