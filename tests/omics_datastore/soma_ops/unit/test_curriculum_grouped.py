"""
Unit tests for curriculum_grouped module.
"""

import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.soma_ops import (
    SomaPrepareCurriculumMetadataError,
    SomaReadError,
    curriculum_grouped,
)
from tests.omics_datastore.soma_ops.conftest import FakeSomaExperiment

# --- Tests for _build_group_filter ---


def test_build_group_filter_single_string_column() -> None:
    """Verify filter for single string column."""
    result = curriculum_grouped._build_group_filter(["tissue"], ("lung",))
    assert result == 'tissue == "lung"'


def test_build_group_filter_multiple_columns() -> None:
    """Verify filter for multiple columns."""
    result = curriculum_grouped._build_group_filter(["tissue", "donor"], ("lung", "D1"))
    assert result == 'tissue == "lung" and donor == "D1"'


def test_build_group_filter_numeric_value() -> None:
    """Verify filter handles numeric values without quotes."""
    result = curriculum_grouped._build_group_filter(["age"], (30,))
    assert result == "age == 30"


def test_build_group_filter_escapes_quotes() -> None:
    """Verify filter escapes quotes in string values."""
    result = curriculum_grouped._build_group_filter(["name"], ('O"Brien',))
    assert result == 'name == "O\\"Brien"'


# --- Tests for prepare_grouped_curriculum validation ---


@pytest.mark.parametrize("bin_size", [0, -10])
def test_prepare_grouped_curriculum_invalid_bin_size(bin_size: int) -> None:
    """Verify non-positive bin_size raises ValueError."""
    with pytest.raises(ValueError, match="bin_size must be positive"):
        curriculum_grouped.prepare_grouped_curriculum(
            experiment_uri="gs://bucket/soma",
            filters=None,
            extract_bin_keys=["tissue"],
            bin_size=bin_size,
        )


def test_prepare_grouped_curriculum_empty_extract_bin_keys() -> None:
    """Verify empty extract_bin_keys raises ValueError."""
    with pytest.raises(ValueError, match="extract_bin_keys cannot be empty"):
        curriculum_grouped.prepare_grouped_curriculum(
            experiment_uri="gs://bucket/soma",
            filters=None,
            extract_bin_keys=[],
            bin_size=100,
        )


def test_prepare_grouped_curriculum_no_cells_raises_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify empty result raises SomaPrepareCurriculumMetadataError."""
    obs_df = pd.DataFrame({"soma_joinid": pd.Series([], dtype="int64"), "tissue": pd.Series([], dtype="str")})

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    with pytest.raises(SomaPrepareCurriculumMetadataError, match="No cells found"):
        curriculum_grouped.prepare_grouped_curriculum(
            experiment_uri="gs://bucket/soma",
            filters=None,
            extract_bin_keys=["tissue"],
            bin_size=100,
        )


def test_prepare_grouped_curriculum_soma_read_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify SOMA read errors are wrapped in SomaReadError."""

    def fake_open(uri: str, mode: str) -> object:
        raise RuntimeError("SOMA connection failed")

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    with pytest.raises(SomaReadError):
        curriculum_grouped.prepare_grouped_curriculum(
            experiment_uri="gs://bucket/soma",
            filters=None,
            extract_bin_keys=["tissue"],
            bin_size=100,
        )


# --- Tests for prepare_grouped_curriculum grouping logic ---


def test_prepare_grouped_curriculum_single_group_single_bin(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify single group smaller than bin_size creates one bin."""
    obs_df = pd.DataFrame({
        "soma_joinid": [1, 2, 3, 4, 5],
        "tissue": ["lung", "lung", "lung", "lung", "lung"],
    })

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = curriculum_grouped.prepare_grouped_curriculum(
        experiment_uri="gs://bucket/soma",
        filters=None,
        extract_bin_keys=["tissue"],
        bin_size=10,
    )

    assert result.num_grouped_bins == 1
    assert result.total_cells == 5
    assert result.grouped_bins[0].group_key == "lung"
    assert result.grouped_bins[0].cell_count == 5
    assert result.grouped_bins[0].joinid_min == 1
    assert result.grouped_bins[0].joinid_max == 5


def test_prepare_grouped_curriculum_single_group_multiple_bins(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify large group is split into multiple bins."""
    obs_df = pd.DataFrame({
        "soma_joinid": [1, 2, 3, 4, 5, 6, 7],
        "tissue": ["lung"] * 7,
    })

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = curriculum_grouped.prepare_grouped_curriculum(
        experiment_uri="gs://bucket/soma",
        filters=None,
        extract_bin_keys=["tissue"],
        bin_size=3,
    )

    assert result.num_grouped_bins == 3
    assert result.total_cells == 7

    # First bin: joinids 1, 2, 3
    assert result.grouped_bins[0].cell_count == 3
    assert result.grouped_bins[0].joinid_min == 1
    assert result.grouped_bins[0].joinid_max == 3

    # Second bin: joinids 4, 5, 6
    assert result.grouped_bins[1].cell_count == 3
    assert result.grouped_bins[1].joinid_min == 4
    assert result.grouped_bins[1].joinid_max == 6

    # Third bin: joinid 7
    assert result.grouped_bins[2].cell_count == 1
    assert result.grouped_bins[2].joinid_min == 7
    assert result.grouped_bins[2].joinid_max == 7


def test_prepare_grouped_curriculum_multiple_groups(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify multiple groups create separate bins."""
    obs_df = pd.DataFrame({
        "soma_joinid": [1, 2, 3, 10, 11, 12],
        "tissue": ["lung", "lung", "lung", "heart", "heart", "heart"],
    })

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = curriculum_grouped.prepare_grouped_curriculum(
        experiment_uri="gs://bucket/soma",
        filters=None,
        extract_bin_keys=["tissue"],
        bin_size=10,
    )

    assert result.num_grouped_bins == 2
    assert result.total_cells == 6

    # Groups are sorted alphabetically: heart, lung
    assert result.grouped_bins[0].group_key == "heart"
    assert result.grouped_bins[0].cell_count == 3

    assert result.grouped_bins[1].group_key == "lung"
    assert result.grouped_bins[1].cell_count == 3


def test_prepare_grouped_curriculum_multiple_grouping_columns(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify grouping by multiple columns creates composite group keys."""
    obs_df = pd.DataFrame({
        "soma_joinid": [1, 2, 3, 4, 5, 6],
        "tissue": ["lung", "lung", "lung", "heart", "heart", "heart"],
        "donor": ["D1", "D1", "D2", "D1", "D1", "D2"],
    })

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = curriculum_grouped.prepare_grouped_curriculum(
        experiment_uri="gs://bucket/soma",
        filters=None,
        extract_bin_keys=["tissue", "donor"],
        bin_size=10,
    )

    # 4 unique combinations: heart/D1, heart/D2, lung/D1, lung/D2
    assert result.num_grouped_bins == 4
    assert result.total_cells == 6

    # Verify group keys are composite
    group_keys = [b.group_key for b in result.grouped_bins]
    assert "heart||D1" in group_keys
    assert "heart||D2" in group_keys
    assert "lung||D1" in group_keys
    assert "lung||D2" in group_keys


def test_prepare_grouped_curriculum_group_filter_format(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify group_filter is correctly formatted for SOMA queries."""
    obs_df = pd.DataFrame({
        "soma_joinid": [1, 2],
        "tissue": ["lung", "heart"],
    })

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = curriculum_grouped.prepare_grouped_curriculum(
        experiment_uri="gs://bucket/soma",
        filters=None,
        extract_bin_keys=["tissue"],
        bin_size=10,
    )

    # Verify group_filter is valid SOMA filter syntax
    assert result.grouped_bins[0].group_filter == 'tissue == "heart"'
    assert result.grouped_bins[1].group_filter == 'tissue == "lung"'


def test_prepare_grouped_curriculum_metadata_fields(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify all metadata fields are correctly populated."""
    obs_df = pd.DataFrame({
        "soma_joinid": [1, 2, 3],
        "tissue": ["lung", "lung", "lung"],
    })

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=obs_df)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    result = curriculum_grouped.prepare_grouped_curriculum(
        experiment_uri="gs://bucket/soma",
        filters={"tissue__eq": "lung"},
        extract_bin_keys=["tissue"],
        bin_size=10,
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        x_layer="raw",
    )

    assert result.experiment_uri == "gs://bucket/soma"
    assert result.filters == {"tissue__eq": "lung"}
    assert result.extract_bin_keys == ["tissue"]
    assert result.obs_columns == ["cell_type"]
    assert result.var_columns == ["symbol"]
    assert result.x_layer == "raw"
    assert result.grouped_bins is not None
    assert result.num_grouped_bins == 1

    # Randomized extraction fields should be None
    assert result.id_ranges is None
    assert result.range_size is None
    assert result.num_ranges is None
