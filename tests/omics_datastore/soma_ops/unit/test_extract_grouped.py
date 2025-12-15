"""
Unit tests for extract_grouped module.
"""

from pathlib import Path

import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.soma_ops import SomaExtractError, extract_grouped
from cellarium.nexus.shared.schemas.omics_datastore import GroupedBin, SomaCurriculumMetadata
from tests.omics_datastore.soma_ops.conftest import FakeAnnData, FakeExecutor, FakeFuture, FakeSomaExperiment


def test_extract_grouped_bin_to_anndata_happy_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify extract_grouped_bin_to_anndata reads SOMA and writes AnnData correctly."""
    experiment_uri = "gs://bucket/soma"
    grouped_bin = GroupedBin(
        group_key="lung",
        group_filter='tissue == "lung"',
        joinid_min=10,
        joinid_max=20,
        cell_count=3,
    )
    output_path = tmp_path / "output.h5ad"

    obs_df = pd.DataFrame({
        "soma_joinid": [10, 15, 20],
        "cell_type": ["T", "B", "T"],
        "tissue": ["lung", "lung", "lung"],
    })

    var_df = pd.DataFrame({
        "soma_joinid": [0, 1, 2],
        "symbol": ["GeneA", "GeneB", "GeneC"],
    })

    x_coo = sp.coo_matrix(
        ([1, 2, 3, 1, 1, 1], ([10, 10, 15, 20, 20, 20], [0, 2, 1, 0, 1, 2])),
        shape=(21, 3),
    )

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        assert uri == experiment_uri
        assert mode == "r"
        return FakeSomaExperiment(obs_data=obs_df, var_data=var_df, x_matrix=x_coo, x_layer="X")

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    adata_instances = []

    def fake_anndata_constructor(X: object = None, obs: object = None, var: object = None) -> FakeAnnData:
        adata = FakeAnnData(X=X, obs=obs, var=var)
        adata_instances.append(adata)
        return adata

    monkeypatch.setattr("anndata.AnnData", fake_anndata_constructor)

    extract_grouped.extract_grouped_bin_to_anndata(
        experiment_uri=experiment_uri,
        grouped_bin=grouped_bin,
        output_path=output_path,
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        x_layer="X",
        output_format="h5ad",
    )

    assert len(adata_instances) == 1
    adata = adata_instances[0]
    assert adata.X.shape == (3, 3)
    assert output_path.exists()


def test_extract_grouped_bin_to_anndata_empty_obs(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify empty obs DataFrame creates empty AnnData file."""
    experiment_uri = "gs://bucket/soma"
    grouped_bin = GroupedBin(
        group_key="mars",
        group_filter='tissue == "mars"',
        joinid_min=10,
        joinid_max=20,
        cell_count=0,
    )
    output_path = tmp_path / "output.h5ad"

    empty_obs = pd.DataFrame({"soma_joinid": pd.Series([], dtype="int64")})

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=empty_obs)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    adata_instances = []

    def fake_anndata_constructor(X: object = None, obs: object = None, var: object = None) -> FakeAnnData:
        adata = FakeAnnData(X=X, obs=obs, var=var)
        adata_instances.append(adata)
        return adata

    monkeypatch.setattr("anndata.AnnData", fake_anndata_constructor)

    extract_grouped.extract_grouped_bin_to_anndata(
        experiment_uri=experiment_uri,
        grouped_bin=grouped_bin,
        output_path=output_path,
        output_format="h5ad",
    )

    assert len(adata_instances) == 1
    assert output_path.exists()


def test_extract_grouped_bin_to_anndata_error_handling(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify SOMA read errors are wrapped in SomaExtractError."""
    experiment_uri = "gs://bucket/soma"
    grouped_bin = GroupedBin(
        group_key="lung",
        group_filter='tissue == "lung"',
        joinid_min=10,
        joinid_max=20,
        cell_count=3,
    )
    output_path = tmp_path / "output.h5ad"

    def fake_open(uri: str, mode: str) -> object:
        raise RuntimeError("SOMA connection failed")

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    with pytest.raises(SomaExtractError):
        extract_grouped.extract_grouped_bin_to_anndata(
            experiment_uri=experiment_uri,
            grouped_bin=grouped_bin,
            output_path=output_path,
            output_format="h5ad",
        )


def test_extract_grouped_bin_worker(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify _extract_grouped_bin_worker calls extract_grouped_bin_to_anndata and returns tuple."""
    extract_calls = []

    def fake_extract(**kwargs: object) -> None:
        extract_calls.append(kwargs)

    monkeypatch.setattr(extract_grouped, "extract_grouped_bin_to_anndata", fake_extract)

    bin_idx = 0
    global_bin_idx = 5
    experiment_uri = "gs://bucket/soma"
    grouped_bin = GroupedBin(
        group_key="lung",
        group_filter='tissue == "lung"',
        joinid_min=10,
        joinid_max=20,
        cell_count=3,
    )
    output_path = tmp_path / "output.h5ad"

    result = extract_grouped._extract_grouped_bin_worker(
        bin_idx=bin_idx,
        global_bin_idx=global_bin_idx,
        experiment_uri=experiment_uri,
        grouped_bin=grouped_bin,
        output_path=output_path,
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        var_joinids=[0, 1, 2],
        x_layer="X",
        output_format="h5ad",
    )

    assert len(extract_calls) == 1
    call_kwargs = extract_calls[0]
    assert call_kwargs["experiment_uri"] == experiment_uri
    assert call_kwargs["grouped_bin"] == grouped_bin
    assert call_kwargs["output_path"] == output_path
    assert call_kwargs["obs_columns"] == ["cell_type"]
    assert call_kwargs["var_columns"] == ["symbol"]
    assert call_kwargs["var_joinids"] == [0, 1, 2]
    assert call_kwargs["x_layer"] == "X"

    assert result == (bin_idx, global_bin_idx, str(output_path))


def test_extract_grouped_bins_no_grouped_bins_raises_error(tmp_path: Path) -> None:
    """Verify extract_grouped_bins raises ValueError when grouped_bins is None."""
    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        total_cells=10,
        grouped_bins=None,
    )

    with pytest.raises(ValueError, match="grouped_bins cannot be None"):
        extract_grouped.extract_grouped_bins(
            curriculum_metadata=curriculum_metadata,
            output_dir=tmp_path,
        )


def test_extract_grouped_bins_happy_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify extract_grouped_bins processes all bins in parallel."""
    grouped_bins = [
        GroupedBin(group_key="lung", group_filter='tissue == "lung"', joinid_min=0, joinid_max=10, cell_count=10),
        GroupedBin(group_key="heart", group_filter='tissue == "heart"', joinid_min=11, joinid_max=20, cell_count=10),
    ]

    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter='tissue in ["lung", "heart"]',
        total_cells=20,
        extract_bin_keys=["tissue"],
        grouped_bins=grouped_bins,
        num_grouped_bins=2,
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        x_layer="X",
    )

    output_dir = tmp_path / "output"

    def fake_worker(
        bin_idx: int,
        global_bin_idx: int,
        experiment_uri: str,
        grouped_bin: GroupedBin,
        output_path: Path,
        obs_columns: list[str] | None,
        var_columns: list[str] | None,
        var_joinids: list[int] | None,
        x_layer: str,
        output_format: str,
    ) -> tuple[int, int, str]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.touch()
        return bin_idx, global_bin_idx, str(output_path)

    monkeypatch.setattr(extract_grouped, "_extract_grouped_bin_worker", fake_worker)

    executor = FakeExecutor(max_workers=2, worker_fn=fake_worker)
    monkeypatch.setattr(extract_grouped, "ProcessPoolExecutor", lambda **kwargs: executor)

    def fake_as_completed(futures: dict[object, int]) -> list[object]:
        return list(futures.keys())

    monkeypatch.setattr(extract_grouped, "as_completed", fake_as_completed)

    extract_grouped.extract_grouped_bins(
        curriculum_metadata=curriculum_metadata,
        output_dir=output_dir,
        max_workers=2,
    )

    assert output_dir.exists()
    assert len(executor.calls) == 2
    assert (output_dir / "extract_000000.h5ad").exists()
    assert (output_dir / "extract_000001.h5ad").exists()


def test_extract_grouped_bins_with_partition_index(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify extract_grouped_bins respects partition_index for distributed execution."""
    grouped_bins = [
        GroupedBin(group_key=f"group_{i}", group_filter=f'id == {i}', joinid_min=i * 10, joinid_max=i * 10 + 9, cell_count=10)
        for i in range(10)
    ]

    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        total_cells=100,
        extract_bin_keys=["id"],
        grouped_bins=grouped_bins,
        num_grouped_bins=10,
    )

    output_dir = tmp_path / "output"

    def fake_worker(
        bin_idx: int,
        global_bin_idx: int,
        experiment_uri: str,
        grouped_bin: GroupedBin,
        output_path: Path,
        obs_columns: list[str] | None,
        var_columns: list[str] | None,
        var_joinids: list[int] | None,
        x_layer: str,
        output_format: str,
    ) -> tuple[int, int, str]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.touch()
        return bin_idx, global_bin_idx, str(output_path)

    monkeypatch.setattr(extract_grouped, "_extract_grouped_bin_worker", fake_worker)

    executor = FakeExecutor(max_workers=2, worker_fn=fake_worker)
    monkeypatch.setattr(extract_grouped, "ProcessPoolExecutor", lambda **kwargs: executor)

    def fake_as_completed(futures: dict[object, int]) -> list[object]:
        return list(futures.keys())

    monkeypatch.setattr(extract_grouped, "as_completed", fake_as_completed)

    # Process partition 1 with 4 bins per partition (bins 4-7)
    extract_grouped.extract_grouped_bins(
        curriculum_metadata=curriculum_metadata,
        output_dir=output_dir,
        partition_index=1,
        max_bins_per_partition=4,
        max_workers=2,
    )

    assert output_dir.exists()
    assert len(executor.calls) == 4

    # Verify output files are named with global bin indices (4, 5, 6, 7)
    assert (output_dir / "extract_000004.h5ad").exists()
    assert (output_dir / "extract_000005.h5ad").exists()
    assert (output_dir / "extract_000006.h5ad").exists()
    assert (output_dir / "extract_000007.h5ad").exists()


def test_extract_grouped_bins_failure_handling(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify extract_grouped_bins raises SomaExtractError when workers fail."""
    grouped_bins = [
        GroupedBin(group_key="lung", group_filter='tissue == "lung"', joinid_min=0, joinid_max=10, cell_count=10),
        GroupedBin(group_key="heart", group_filter='tissue == "heart"', joinid_min=11, joinid_max=20, cell_count=10),
    ]

    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        total_cells=20,
        extract_bin_keys=["tissue"],
        grouped_bins=grouped_bins,
        num_grouped_bins=2,
    )

    output_dir = tmp_path / "output"

    class FailingExecutor(FakeExecutor):
        def submit(self, fn: object, *args: object, **kwargs: object) -> FakeFuture:
            bin_idx = args[0]
            if bin_idx == 1:
                future = FakeFuture(exception=RuntimeError(f"Worker {bin_idx} failed"))
            else:
                future = FakeFuture(result=(bin_idx, bin_idx, f"path_{bin_idx}"))
            self.futures.append(future)
            return future

    executor = FailingExecutor(max_workers=2)
    monkeypatch.setattr(extract_grouped, "ProcessPoolExecutor", lambda **kwargs: executor)

    def fake_as_completed(futures: dict[object, int]) -> list[object]:
        return list(futures.keys())

    monkeypatch.setattr(extract_grouped, "as_completed", fake_as_completed)

    with pytest.raises(SomaExtractError):
        extract_grouped.extract_grouped_bins(
            curriculum_metadata=curriculum_metadata,
            output_dir=output_dir,
        )


def test_extract_grouped_bins_empty_partition(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Verify extract_grouped_bins handles empty partition gracefully."""
    grouped_bins = [
        GroupedBin(group_key="lung", group_filter='tissue == "lung"', joinid_min=0, joinid_max=10, cell_count=10),
    ]

    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        total_cells=10,
        extract_bin_keys=["tissue"],
        grouped_bins=grouped_bins,
        num_grouped_bins=1,
    )

    output_dir = tmp_path / "output"

    # Partition 5 with 2 bins per partition -> no bins to process (only 1 bin total)
    extract_grouped.extract_grouped_bins(
        curriculum_metadata=curriculum_metadata,
        output_dir=output_dir,
        partition_index=5,
        max_bins_per_partition=2,
    )

    # Should complete without error, no files created
    assert not output_dir.exists() or len(list(output_dir.iterdir())) == 0
