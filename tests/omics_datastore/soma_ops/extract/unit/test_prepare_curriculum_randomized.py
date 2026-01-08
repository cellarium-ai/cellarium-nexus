import numpy as np
import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.soma_ops import (
    SomaFilterError,
    SomaPrepareCurriculumMetadataError,
    SomaReadError,
)
from cellarium.nexus.omics_datastore.soma_ops._extract import prepare_curriculum_randomized
from cellarium.nexus.shared.schemas.omics_datastore import IdContiguousRange
from tests.omics_datastore.soma_ops.conftest import FakeSomaExperiment


def test_compute_contiguous_ranges_basic_chunking() -> None:
    """
    Verify basic chunking of sorted values into contiguous ranges.
    """
    values = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    ranges = prepare_curriculum_randomized.compute_contiguous_ranges(values=values, chunk_size=2, shuffle=False)
    assert ranges == [
        IdContiguousRange(start=1, end=2),
        IdContiguousRange(start=3, end=4),
        IdContiguousRange(start=5, end=5),
    ], "Expected 3 ranges with correct start/end values"


def test_compute_contiguous_ranges_non_consecutive_values() -> None:
    """
    Verify that non-consecutive values are still chunked by position.
    """
    values = np.array([10, 20, 21, 100], dtype=np.int64)
    ranges = prepare_curriculum_randomized.compute_contiguous_ranges(values=values, chunk_size=2, shuffle=False)
    assert ranges == [
        IdContiguousRange(start=10, end=20),
        IdContiguousRange(start=21, end=100),
    ], "Expected 2 ranges chunked by position"


def test_compute_contiguous_ranges_empty_values() -> None:
    """
    Verify empty values array returns empty ranges list.
    """
    values = np.array([], dtype=np.int64)
    ranges = prepare_curriculum_randomized.compute_contiguous_ranges(values=values, chunk_size=10, shuffle=False)
    assert ranges == [], "Expected empty ranges list for empty input"


def test_shuffle_but_keep_last():
    """
    Test `shuffle_but_keep_last` function
    """
    values = [1, 2, 3, 4, 99]
    values_original = list(values)

    out = prepare_curriculum_randomized.shuffle_but_keep_last(values=values)

    assert out[-1] == values[-1], "Last element must stay the same"
    assert len(out) == len(values), "Lengths must stay the same"
    assert sorted(out[:-1]) == sorted(values[:-1]), "Prefix elements must be a permutation of the original prefix"
    assert values_original == values, "Function must be immutable: original values must stay the same"


@pytest.mark.parametrize("chunk_size", [0, -1])
def test_compute_contiguous_ranges_invalid_chunk_size(chunk_size: int) -> None:
    """
    Verify non-positive chunk_size raises ValueError.
    """
    values = np.array([1, 2, 3], dtype=np.int64)
    with pytest.raises(ValueError):
        prepare_curriculum_randomized.compute_contiguous_ranges(values=values, chunk_size=chunk_size)


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

    result = prepare_curriculum_randomized.read_filtered_joinids(
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

    result = prepare_curriculum_randomized.read_filtered_joinids(
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
        prepare_curriculum_randomized.read_filtered_joinids(
            experiment_uri=experiment_uri,
            value_filter=value_filter,
        )


@pytest.mark.parametrize("range_size", [0, -5])
def test_plan_soma_extract_invalid_range_size(range_size: int) -> None:
    """
    Verify non-positive range_size raises ValueError.
    """
    with pytest.raises(ValueError):
        prepare_curriculum_randomized.prepare_extract_curriculum(
            experiment_uri="gs://bucket/soma",
            filters=None,
            range_size=range_size,
            extract_bin_size=100,
        )


def test_plan_soma_extract_no_cells(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify plan with no matching cells returns empty ranges.
    """
    experiment_uri = "gs://bucket/soma"
    filters = {"organism__eq": "Martian"}

    # Mock read_filtered_joinids to return empty array
    monkeypatch.setattr(
        prepare_curriculum_randomized,
        "read_filtered_joinids",
        lambda experiment_uri, value_filter: np.array([], dtype=np.int64),
    )
    with pytest.raises(SomaPrepareCurriculumMetadataError):
        _ = prepare_curriculum_randomized.prepare_extract_curriculum(
            experiment_uri=experiment_uri,
            filters=filters,
            range_size=100,
            extract_bin_size=100,
            shuffle_ranges=False,
        )


def test_plan_soma_extract_normal_flow(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify plan_soma_extract creates correct ranges without shuffling.
    """
    experiment_uri = "gs://bucket/soma"
    filters = {"tissue__eq": "lung"}

    # Mock read_filtered_joinids
    monkeypatch.setattr(
        prepare_curriculum_randomized,
        "read_filtered_joinids",
        lambda experiment_uri, value_filter: np.array([5, 10, 15, 20], dtype=np.int64),
    )

    plan = prepare_curriculum_randomized.prepare_extract_curriculum(
        experiment_uri=experiment_uri,
        filters=filters,
        range_size=2,
        extract_bin_size=2,
        shuffle_ranges=False,
    )

    assert plan.experiment_uri == experiment_uri, "experiment_uri should match input"
    assert plan.value_filter == '(tissue == "lung")', "value_filter should be translated correctly"
    assert plan.total_cells == 4, "Expected 4 total cells"
    assert plan.range_size == 2, "range_size should match input"
    assert plan.extract_bin_size == 2, "extract_bin_size should match input"
    assert plan.filters == filters, "filters should match input"
    assert len(plan.id_ranges) == 2, "Expected 2 id_ranges"
    assert plan.id_ranges[0] == IdContiguousRange(start=5, end=10), "First range should be (5, 10)"
    assert plan.id_ranges[1] == IdContiguousRange(start=15, end=20), "Second range should be (15, 20)"
    assert plan.extract_bin_indexes == [0, 1], "Extract bin indexes should be sequential when not shuffled"


def test_plan_soma_extract_with_shuffle(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify shuffle_ranges=True shuffles the ranges list.
    """
    experiment_uri = "gs://bucket/soma"
    filters = None

    # Mock read_filtered_joinids
    monkeypatch.setattr(
        prepare_curriculum_randomized,
        "read_filtered_joinids",
        lambda experiment_uri, value_filter: np.array([1, 2, 3, 4, 5, 6], dtype=np.int64),
    )

    # Mock random.sample to verify it's called
    shuffle_called = []

    def _fake_sample(lst: list[object], k: int) -> list[object]:
        shuffle_called.append(True)
        # Reverse for deterministic test
        return list(reversed(lst))

    import random

    monkeypatch.setattr(random, "sample", _fake_sample)

    plan = prepare_curriculum_randomized.prepare_extract_curriculum(
        experiment_uri=experiment_uri,
        filters=filters,
        range_size=2,
        extract_bin_size=2,
        shuffle_ranges=True,
        shuffle_extract_bin_indexes=True,
    )

    assert (
        len(shuffle_called) == 2
    ), "`random.sample` should be called 2 times (one for ranges, one for shuffling extract bin indexes"
    assert plan.total_cells == 6, "Expected 6 total cells"

    assert len(plan.id_ranges) == 3, "Expected 3 id_ranges"
    assert plan.id_ranges[0] == IdContiguousRange(start=3, end=4), "First range should be (3, 4) after shuffle"
    assert plan.id_ranges[1] == IdContiguousRange(start=1, end=2), "Second range should be (1, 2) after shuffle"
    assert plan.id_ranges[2] == IdContiguousRange(start=5, end=6), "Last range should stay (5, 6)"

    assert plan.extract_bin_indexes != [
        0,
        1,
        2,
    ], "Extract bin indexes should not be sequential when shuffle_extract_bin_indexes=True"


def test_plan_soma_extract_filter_error_propagation(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify SomaFilterError from build_soma_value_filter propagates.
    """
    # Use a filter with an unsupported operator to trigger SomaFilterError
    with pytest.raises(SomaFilterError):
        prepare_curriculum_randomized.prepare_extract_curriculum(
            experiment_uri="gs://bucket/soma",
            filters={"bad__unsupported_op": "value"},
            range_size=100,
            extract_bin_size=100,
        )


def test_plan_soma_extract_read_error_propagation(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify SomaReadError from read_filtered_joinids propagates.
    """
    from cellarium.nexus.omics_datastore.soma_ops import filters as filters_module

    monkeypatch.setattr(filters_module, "build_soma_value_filter", lambda filters: "")

    def _raise_read_error(experiment_uri: str, value_filter: str) -> np.ndarray:
        raise SomaReadError("Read failed")

    monkeypatch.setattr(prepare_curriculum_randomized, "read_filtered_joinids", _raise_read_error)

    with pytest.raises(SomaReadError):
        prepare_curriculum_randomized.prepare_extract_curriculum(
            experiment_uri="gs://bucket/soma",
            filters=None,
            range_size=100,
            extract_bin_size=100,
        )
