from pathlib import Path

import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.soma_ops import SomaReadError
from cellarium.nexus.omics_datastore.soma_ops import data_operator as data_operator_module
from cellarium.nexus.shared.schemas.omics_datastore import IdContiguousRange, RandomizedCurriculumMetadata


def test_init_valid_experiment_uri() -> None:
    """
    Verify TileDBSOMADataOperator initializes with valid experiment_uri.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma_experiment")
    assert operator.experiment_uri == "gs://bucket/soma_experiment"


def test_init_empty_experiment_uri() -> None:
    """
    Verify empty experiment_uri raises ValueError.
    """
    with pytest.raises(ValueError):
        data_operator_module.TileDBSOMADataOperator(experiment_uri="")


def test_count_cells_happy_path(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify count_cells reads obs and returns correct count.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock tiledbsoma.open
    class _FakePandasTable:
        def to_pandas(self) -> pd.DataFrame:
            return pd.DataFrame({"soma_joinid": [1, 2, 3, 4, 5]})

    class _FakeQuery:
        def concat(self) -> _FakePandasTable:
            return _FakePandasTable()

    class _FakeObs:
        def read(self, column_names: list[str], value_filter: str | None) -> _FakeQuery:
            assert column_names == ["soma_joinid"]
            # Filter will be built by real build_soma_value_filter
            assert value_filter == '(tissue == "lung")' or value_filter is None
            return _FakeQuery()

    class _FakeExp:
        def __init__(self) -> None:
            self.obs = _FakeObs()

        def __enter__(self) -> "_FakeExp":
            return self

        def __exit__(self, *args: object) -> None:
            pass

    def _fake_open(uri: str, mode: str) -> _FakeExp:
        assert uri == "gs://bucket/soma"
        assert mode == "r"
        return _FakeExp()

    import tiledbsoma

    monkeypatch.setattr(tiledbsoma, "open", _fake_open)

    # Execute
    filters = {"tissue__eq": "lung"}
    count = operator.count_cells(filter_statements=filters)

    # Verify
    assert count == 5


def test_count_cells_empty_filter(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify count_cells with None filter passes None to SOMA.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock build_soma_value_filter
    from cellarium.nexus.omics_datastore.soma_ops import filters as filters_module

    monkeypatch.setattr(filters_module, "build_soma_value_filter", lambda filters: "")

    # Mock tiledbsoma.open
    class _FakePandasTable:
        def to_pandas(self) -> pd.DataFrame:
            return pd.DataFrame({"soma_joinid": [1, 2, 3]})

    class _FakeQuery:
        def concat(self) -> _FakePandasTable:
            return _FakePandasTable()

    class _FakeObs:
        def read(self, column_names: list[str], value_filter: str | None) -> _FakeQuery:
            assert column_names == ["soma_joinid"]
            assert value_filter is None or value_filter == ""
            return _FakeQuery()

    class _FakeExp:
        def __init__(self) -> None:
            self.obs = _FakeObs()

        def __enter__(self) -> "_FakeExp":
            return self

        def __exit__(self, *args: object) -> None:
            pass

    def _fake_open(uri: str, mode: str) -> _FakeExp:
        return _FakeExp()

    import tiledbsoma

    monkeypatch.setattr(tiledbsoma, "open", _fake_open)

    # Execute
    count = operator.count_cells(filter_statements=None)

    # Verify
    assert count == 3


def test_count_cells_error_handling(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify count_cells wraps SOMA errors in SomaReadError.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock build_soma_value_filter
    from cellarium.nexus.omics_datastore.soma_ops import filters as filters_module

    monkeypatch.setattr(filters_module, "build_soma_value_filter", lambda filters: "")

    # Mock tiledbsoma.open to raise
    def _fake_open(uri: str, mode: str) -> object:
        raise RuntimeError("SOMA connection failed")

    import tiledbsoma

    monkeypatch.setattr(tiledbsoma, "open", _fake_open)

    # Execute and expect error
    with pytest.raises(SomaReadError):
        operator.count_cells(filter_statements=None)


def test_compute_extract_plan_delegates_to_planning(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify compute_extract_plan delegates to plan_soma_extract.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock plan_soma_extract at the data_operator module level
    plan_calls = []

    def _fake_plan_soma_extract(
        experiment_uri: str,
        filters: object,
        range_size: int,
        extract_bin_size: int,
        shuffle_ranges: bool,
        var_filter_column: str | None,
        var_filter_values: list[str] | None,
        obs_columns: list[str] | None,
        var_columns: list[str] | None,
        x_layer: str,
    ) -> RandomizedCurriculumMetadata:
        plan_calls.append(
            {
                "experiment_uri": experiment_uri,
                "filters": filters,
                "range_size": range_size,
                "extract_bin_size": extract_bin_size,
                "shuffle_ranges": shuffle_ranges,
                "var_filter_column": var_filter_column,
                "var_filter_values": var_filter_values,
                "obs_columns": obs_columns,
                "var_columns": var_columns,
                "x_layer": x_layer,
            }
        )
        return RandomizedCurriculumMetadata(
            experiment_uri=experiment_uri,
            value_filter='tissue == "lung"',
            id_ranges=[IdContiguousRange(start=0, end=10)],
            total_cells=10,
            range_size=range_size,
            num_ranges=1,
            extract_bin_size=extract_bin_size,
            num_bins=1,
            last_bin_size=10,
            extract_bin_indexes=[0],
            filters=filters,
            var_filter_column=var_filter_column,
            var_filter_values=var_filter_values,
            obs_columns=obs_columns,
            var_columns=var_columns,
            x_layer=x_layer,
        )

    monkeypatch.setattr(data_operator_module, "prepare_extract_curriculum", _fake_plan_soma_extract)

    # Execute
    filters = {"tissue__eq": "lung"}
    plan = operator.prepare_curriculum_metadata(
        filters=filters,
        range_size=100,
        extract_bin_size=100,
        shuffle_ranges=False,
        var_filter_column="gene_symbol",
        var_filter_values=["ACTB", "GAPDH"],
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        x_layer="raw",
    )

    # Verify delegation
    assert len(plan_calls) == 1
    call = plan_calls[0]
    assert call["experiment_uri"] == "gs://bucket/soma"
    assert call["filters"] == filters
    assert call["range_size"] == 100
    assert call["shuffle_ranges"] is False
    assert call["var_filter_column"] == "gene_symbol"
    assert call["var_filter_values"] == ["ACTB", "GAPDH"]
    assert call["obs_columns"] == ["cell_type"]
    assert call["var_columns"] == ["symbol"]
    assert call["x_layer"] == "raw"

    # Verify return value
    assert plan.experiment_uri == "gs://bucket/soma"
    assert plan.total_cells == 10
    assert plan.range_size == 100
    assert plan.var_filter_column == "gene_symbol"
    assert plan.var_filter_values == ["ACTB", "GAPDH"]
    assert plan.obs_columns == ["cell_type"]
    assert plan.var_columns == ["symbol"]
    assert plan.x_layer == "raw"


def test_extract_randomized_with_temp_dir(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify extract_randomized uses provided temp_dir.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock extract_ranges and shuffle_extracted_chunks at data_operator module level
    extract_calls = []
    shuffle_calls = []

    def _fake_extract_ranges(**kwargs: object) -> None:
        extract_calls.append(kwargs)

    def _fake_shuffle(**kwargs: object) -> None:
        shuffle_calls.append(kwargs)

    monkeypatch.setattr(data_operator_module, "extract_ranges", _fake_extract_ranges)
    monkeypatch.setattr(data_operator_module, "shuffle_extracted_chunks", _fake_shuffle)

    # Mock shutil.rmtree
    import shutil

    rmtree_calls = []

    def _fake_rmtree(path: Path) -> None:
        rmtree_calls.append(path)

    monkeypatch.setattr(shutil, "rmtree", _fake_rmtree)

    # Execute
    plan = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        id_ranges=[IdContiguousRange(start=0, end=10)],
        total_cells=10,
        range_size=10,
        num_ranges=1,
        extract_bin_size=10,
        num_bins=1,
        last_bin_size=10,
        extract_bin_indexes=[0],
        filters=None,
    )

    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()
    output_dir = tmp_path / "output"

    operator.extract_randomized(
        curriculum_metadata=plan,
        output_dir=output_dir,
        output_format="zarr",
        temp_dir=temp_dir,
        max_workers_extract=2,
        max_workers_shuffle=4,
        cleanup_temp=True,
    )

    # Verify extract_ranges was called with temp_dir
    assert len(extract_calls) == 1
    assert extract_calls[0]["output_dir"] == temp_dir
    assert extract_calls[0]["curriculum_metadata"] == plan
    assert extract_calls[0]["partition_index"] == 0
    assert extract_calls[0]["max_ranges_per_partition"] is None
    assert extract_calls[0]["max_workers"] == 2

    # Verify shuffle was called
    assert len(shuffle_calls) == 1
    assert shuffle_calls[0]["input_dir"] == temp_dir
    assert shuffle_calls[0]["output_dir"] == output_dir
    assert shuffle_calls[0]["curriculum_metadata"] == plan
    assert shuffle_calls[0]["partition_index"] == 0
    assert shuffle_calls[0]["max_output_chunks_per_partition"] is None
    assert shuffle_calls[0]["max_workers"] == 4

    # Verify cleanup was called
    assert rmtree_calls == [temp_dir]


def test_extract_randomized_auto_temp_dir(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify extract_randomized creates temp_dir when None.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock extract_ranges and shuffle_extracted_chunks at data_operator module level
    extract_calls = []
    shuffle_calls = []

    def _fake_extract_ranges(**kwargs: object) -> None:
        extract_calls.append(kwargs)

    def _fake_shuffle(**kwargs: object) -> None:
        shuffle_calls.append(kwargs)

    monkeypatch.setattr(data_operator_module, "extract_ranges", _fake_extract_ranges)
    monkeypatch.setattr(data_operator_module, "shuffle_extracted_chunks", _fake_shuffle)

    # Mock tempfile.mkdtemp
    import tempfile

    mkdtemp_calls = []

    def _fake_mkdtemp(prefix: str) -> str:
        mkdtemp_calls.append(prefix)
        temp = tmp_path / "auto_temp"
        temp.mkdir()
        return str(temp)

    monkeypatch.setattr(tempfile, "mkdtemp", _fake_mkdtemp)

    # Mock shutil.rmtree
    import shutil

    rmtree_calls = []

    def _fake_rmtree(path: Path) -> None:
        rmtree_calls.append(path)

    monkeypatch.setattr(shutil, "rmtree", _fake_rmtree)

    # Execute
    plan = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        id_ranges=[IdContiguousRange(start=0, end=10)],
        total_cells=10,
        range_size=10,
        num_ranges=1,
        extract_bin_size=10,
        num_bins=1,
        last_bin_size=10,
        extract_bin_indexes=[0],
        filters=None,
    )

    output_dir = tmp_path / "output"

    operator.extract_randomized(
        curriculum_metadata=plan,
        output_dir=output_dir,
        output_format="h5ad",
        temp_dir=None,
        cleanup_temp=True,
    )

    # Verify mkdtemp was called
    assert len(mkdtemp_calls) == 1
    assert mkdtemp_calls[0] == "soma_extract_temp_"

    # Verify extract_ranges was called with auto temp_dir
    assert len(extract_calls) == 1
    assert str(extract_calls[0]["output_dir"]) == str(tmp_path / "auto_temp")

    # Verify cleanup was called
    assert len(rmtree_calls) == 1


def test_extract_randomized_no_cleanup(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify extract_randomized skips cleanup when cleanup_temp=False.
    """
    operator = data_operator_module.TileDBSOMADataOperator(experiment_uri="gs://bucket/soma")

    # Mock extract_ranges and shuffle_extracted_chunks at data_operator module level
    def _fake_extract_ranges(**kwargs: object) -> None:
        pass

    def _fake_shuffle(**kwargs: object) -> None:
        pass

    monkeypatch.setattr(data_operator_module, "extract_ranges", _fake_extract_ranges)
    monkeypatch.setattr(data_operator_module, "shuffle_extracted_chunks", _fake_shuffle)

    # Mock shutil.rmtree
    import shutil

    rmtree_calls = []

    def _fake_rmtree(path: Path) -> None:
        rmtree_calls.append(path)

    monkeypatch.setattr(shutil, "rmtree", _fake_rmtree)

    # Execute
    plan = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        id_ranges=[IdContiguousRange(start=0, end=10)],
        total_cells=10,
        range_size=10,
        num_ranges=1,
        extract_bin_size=10,
        num_bins=1,
        last_bin_size=10,
        extract_bin_indexes=[0],
        filters=None,
    )

    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()
    output_dir = tmp_path / "output"

    operator.extract_randomized(
        curriculum_metadata=plan,
        output_dir=output_dir,
        temp_dir=temp_dir,
        cleanup_temp=False,
    )

    # Verify cleanup was NOT called
    assert rmtree_calls == []
