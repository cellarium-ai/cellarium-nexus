import numpy as np
import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.soma_ops import SomaFilterError, SomaReadError
from cellarium.nexus.omics_datastore.soma_ops import planning as planning_module
from cellarium.nexus.shared.schemas.omics_datastore import SomaJoinIdRange
from tests.omics_datastore.soma_ops.conftest import FakeSomaExperiment


def test_compute_contiguous_ranges_basic_chunking() -> None:
    """
    Verify basic chunking of sorted values into contiguous ranges.
    """
    values = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    ranges = planning_module.compute_contiguous_ranges(values=values, chunk_size=2)
    assert ranges == [
        SomaJoinIdRange(start=1, end=2),
        SomaJoinIdRange(start=3, end=4),
        SomaJoinIdRange(start=5, end=5),
    ]


def test_compute_contiguous_ranges_non_consecutive_values() -> None:
    """
    Verify that non-consecutive values are still chunked by position.
    """
    values = np.array([10, 20, 21, 100], dtype=np.int64)
    ranges = planning_module.compute_contiguous_ranges(values=values, chunk_size=2)
    assert ranges == [
        SomaJoinIdRange(start=10, end=20),
        SomaJoinIdRange(start=21, end=100),
    ]


def test_compute_contiguous_ranges_empty_values() -> None:
    """
    Verify empty values array returns empty ranges list.
    """
    values = np.array([], dtype=np.int64)
    ranges = planning_module.compute_contiguous_ranges(values=values, chunk_size=10)
    assert ranges == []


@pytest.mark.parametrize("chunk_size", [0, -1])
def test_compute_contiguous_ranges_invalid_chunk_size(chunk_size: int) -> None:
    """
    Verify non-positive chunk_size raises ValueError.
    """
    values = np.array([1, 2, 3], dtype=np.int64)
    with pytest.raises(ValueError):
        planning_module.compute_contiguous_ranges(values=values, chunk_size=chunk_size)


def test_read_filtered_joinids_happy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify read_filtered_joinids reads and sorts soma_joinid column.
    """
    experiment_uri = "gs://bucket/soma_experiment"
    value_filter = 'tissue == "lung"'

    # Return unsorted joinids
    obs_data = pd.DataFrame({"soma_joinid": [30, 10, 20, 5]})

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        assert uri == experiment_uri
        assert mode == "r"
        return FakeSomaExperiment(obs_data, expected_filter='tissue == "lung"')

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = planning_module.read_filtered_joinids(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
    )

    # Should be sorted
    expected = np.array([5, 10, 20, 30], dtype=np.int64)
    np.testing.assert_array_equal(result, expected)


def test_read_filtered_joinids_empty_filter(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify empty value_filter is passed as None to SOMA.
    """
    experiment_uri = "gs://bucket/soma_experiment"
    value_filter = ""

    obs_data = pd.DataFrame({"soma_joinid": [1, 2, 3]})

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = planning_module.read_filtered_joinids(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
    )

    expected = np.array([1, 2, 3], dtype=np.int64)
    np.testing.assert_array_equal(result, expected)


def test_read_filtered_joinids_error_handling(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify SOMA read errors are wrapped in SomaReadError.
    """
    experiment_uri = "gs://bucket/soma_experiment"
    value_filter = 'tissue == "lung"'

    def fake_open(uri: str, mode: str) -> object:
        raise RuntimeError("SOMA connection failed")

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    with pytest.raises(SomaReadError):
        planning_module.read_filtered_joinids(
            experiment_uri=experiment_uri,
            value_filter=value_filter,
        )


@pytest.mark.parametrize("range_size", [0, -5])
def test_plan_soma_extract_invalid_range_size(range_size: int) -> None:
    """
    Verify non-positive range_size raises ValueError.
    """
    with pytest.raises(ValueError):
        planning_module.plan_soma_extract(
            experiment_uri="gs://bucket/soma",
            filters=None,
            range_size=range_size,
        )


def test_plan_soma_extract_no_cells(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify plan with no matching cells returns empty ranges.
    """
    experiment_uri = "gs://bucket/soma"
    filters = {"organism__eq": "Martian"}

    # Mock read_filtered_joinids to return empty array
    monkeypatch.setattr(
        planning_module,
        "read_filtered_joinids",
        lambda experiment_uri, value_filter: np.array([], dtype=np.int64),
    )

    plan = planning_module.plan_soma_extract(
        experiment_uri=experiment_uri,
        filters=filters,
        range_size=100,
        shuffle_ranges=False,
    )

    assert plan.experiment_uri == experiment_uri
    assert plan.value_filter == '(organism == "Martian")'
    assert plan.joinid_ranges == []
    assert plan.total_cells == 0
    assert plan.range_size == 100
    assert plan.filters == filters


def test_plan_soma_extract_normal_flow(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify plan_soma_extract creates correct ranges without shuffling.
    """
    experiment_uri = "gs://bucket/soma"
    filters = {"tissue__eq": "lung"}

    # Mock read_filtered_joinids
    monkeypatch.setattr(
        planning_module,
        "read_filtered_joinids",
        lambda experiment_uri, value_filter: np.array([5, 10, 15, 20], dtype=np.int64),
    )

    plan = planning_module.plan_soma_extract(
        experiment_uri=experiment_uri,
        filters=filters,
        range_size=2,
        shuffle_ranges=False,
    )

    assert plan.experiment_uri == experiment_uri
    assert plan.value_filter == '(tissue == "lung")'
    assert plan.total_cells == 4
    assert plan.range_size == 2
    assert plan.filters == filters
    assert len(plan.joinid_ranges) == 2
    assert plan.joinid_ranges[0] == SomaJoinIdRange(start=5, end=10)
    assert plan.joinid_ranges[1] == SomaJoinIdRange(start=15, end=20)


def test_plan_soma_extract_with_shuffle(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify shuffle_ranges=True shuffles the ranges list.
    """
    experiment_uri = "gs://bucket/soma"
    filters = None

    # Mock read_filtered_joinids
    monkeypatch.setattr(
        planning_module,
        "read_filtered_joinids",
        lambda experiment_uri, value_filter: np.array([1, 2, 3, 4, 5, 6], dtype=np.int64),
    )

    # Mock random.shuffle to verify it's called
    shuffle_called = []

    def _fake_shuffle(lst: list[object]) -> None:
        shuffle_called.append(True)
        # Reverse for deterministic test
        lst.reverse()

    import random

    monkeypatch.setattr(random, "shuffle", _fake_shuffle)

    plan = planning_module.plan_soma_extract(
        experiment_uri=experiment_uri,
        filters=filters,
        range_size=2,
        shuffle_ranges=True,
    )

    assert len(shuffle_called) == 1
    assert plan.total_cells == 6
    # Ranges should be reversed due to our fake shuffle
    assert len(plan.joinid_ranges) == 3
    # After reverse: last range becomes first
    assert plan.joinid_ranges[0] == SomaJoinIdRange(start=5, end=6)


def test_plan_soma_extract_filter_error_propagation(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify SomaFilterError from build_soma_value_filter propagates.
    """
    # Use a filter with an unsupported operator to trigger SomaFilterError
    with pytest.raises(SomaFilterError):
        planning_module.plan_soma_extract(
            experiment_uri="gs://bucket/soma",
            filters={"bad__unsupported_op": "value"},
            range_size=100,
        )


def test_plan_soma_extract_read_error_propagation(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify SomaReadError from read_filtered_joinids propagates.
    """
    from cellarium.nexus.omics_datastore.soma_ops import filters as filters_module

    monkeypatch.setattr(filters_module, "build_soma_value_filter", lambda filters: "")

    def _raise_read_error(experiment_uri: str, value_filter: str) -> np.ndarray:
        raise SomaReadError("Read failed")

    monkeypatch.setattr(planning_module, "read_filtered_joinids", _raise_read_error)

    with pytest.raises(SomaReadError):
        planning_module.plan_soma_extract(
            experiment_uri="gs://bucket/soma",
            filters=None,
            range_size=100,
        )
