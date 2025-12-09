from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.soma_ops import SomaExtractError, extract
from cellarium.nexus.shared.schemas.omics_datastore import IdContiguousRange, SomaCurriculumMetadata
from tests.omics_datastore.soma_ops.conftest import FakeAnnData, FakeExecutor, FakeFuture, FakeSomaExperiment


def test_extract_range_to_anndata_happy_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify extract_range_to_anndata reads SOMA and writes AnnData correctly.
    """
    experiment_uri = "gs://bucket/soma"
    value_filter = 'tissue == "lung"'
    joinid_range = IdContiguousRange(start=10, end=20)
    output_path = tmp_path / "output.h5ad"

    # Mock SOMA experiment
    obs_df = pd.DataFrame(
        {
            "soma_joinid": [10, 15, 20],
            "cell_type": ["T", "B", "T"],
            "tissue": ["lung", "lung", "lung"],
        }
    )

    var_df = pd.DataFrame(
        {
            "soma_joinid": [0, 1, 2],
            "symbol": ["GeneA", "GeneB", "GeneC"],
        }
    )

    # Create sparse matrix with soma_joinids as row indices
    # Row indices should match soma_joinids: [10, 15, 20]
    x_coo = sp.coo_matrix(
        ([1, 2, 3, 1, 1, 1], ([10, 10, 15, 20, 20, 20], [0, 2, 1, 0, 1, 2])),
        shape=(21, 3),  # Max joinid is 20, so shape is 21 x 3
    )

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        assert uri == experiment_uri
        assert mode == "r"
        return FakeSomaExperiment(obs_data=obs_df, var_data=var_df, x_matrix=x_coo, x_layer="X")

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    # Mock AnnData
    adata_instances = []

    def fake_anndata_constructor(X: object = None, obs: object = None, var: object = None) -> FakeAnnData:
        adata = FakeAnnData(X=X, obs=obs, var=var)
        adata_instances.append(adata)
        return adata

    monkeypatch.setattr("anndata.AnnData", fake_anndata_constructor)

    # Execute
    extract.extract_range_to_anndata(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        joinid_range=joinid_range,
        output_path=output_path,
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        x_layer="X",
        output_format="h5ad",
    )

    # Verify AnnData was created
    assert len(adata_instances) == 1
    adata = adata_instances[0]
    # Verify X matrix shape is correct (remapped to 0-based indices)
    assert adata.X.shape == (3, 3)
    assert adata.X.nnz == 6
    # Verify obs index is soma_joinid
    assert list(adata.obs.index) == [10, 15, 20]
    # Verify var index is soma_joinid
    assert list(adata.var.index) == [0, 1, 2]
    # Verify output file was created
    assert output_path.exists()


def test_extract_range_to_anndata_empty_obs(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify empty obs DataFrame creates empty AnnData file.
    """
    experiment_uri = "gs://bucket/soma"
    value_filter = 'tissue == "mars"'
    joinid_range = IdContiguousRange(start=10, end=20)
    output_path = tmp_path / "output.h5ad"

    # Mock SOMA experiment returning empty obs
    empty_obs = pd.DataFrame({"soma_joinid": []})

    def fake_open(uri: str, mode: str) -> FakeSomaExperiment:
        return FakeSomaExperiment(obs_data=empty_obs)

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    # Mock AnnData
    adata_instances = []

    def fake_anndata_constructor(X: object = None, obs: object = None, var: object = None) -> FakeAnnData:
        adata = FakeAnnData(X=X, obs=obs, var=var)
        adata_instances.append(adata)
        return adata

    monkeypatch.setattr("anndata.AnnData", fake_anndata_constructor)

    # Execute
    extract.extract_range_to_anndata(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        joinid_range=joinid_range,
        output_path=output_path,
        output_format="h5ad",
    )

    # Verify empty AnnData was created
    assert len(adata_instances) == 1
    assert output_path.exists()


def test_extract_range_to_anndata_error_handling(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify SOMA read errors are wrapped in SomaExtractError.
    """
    experiment_uri = "gs://bucket/soma"
    value_filter = ""
    joinid_range = IdContiguousRange(start=10, end=20)
    output_path = tmp_path / "output.h5ad"

    def fake_open(uri: str, mode: str) -> object:
        raise RuntimeError("SOMA connection failed")

    monkeypatch.setattr("tiledbsoma.open", fake_open)

    with pytest.raises(SomaExtractError):
        extract.extract_range_to_anndata(
            experiment_uri=experiment_uri,
            value_filter=value_filter,
            joinid_range=joinid_range,
            output_path=output_path,
            output_format="h5ad",
        )


def test_extract_range_worker(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify _extract_range_worker calls extract_range_to_anndata and returns tuple.
    """
    extract_calls = []

    def fake_extract(**kwargs: object) -> None:
        extract_calls.append(kwargs)

    monkeypatch.setattr(extract, "extract_range_to_anndata", fake_extract)

    idx = 5
    experiment_uri = "gs://bucket/soma"
    value_filter = 'tissue == "lung"'
    joinid_range = IdContiguousRange(start=10, end=20)
    output_path = tmp_path / "output.h5ad"
    obs_columns = ["cell_type"]
    var_columns = ["symbol"]
    x_layer = "raw"

    result = extract._extract_range_worker(
        idx=idx,
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        joinid_range=joinid_range,
        output_path=output_path,
        obs_columns=obs_columns,
        var_columns=var_columns,
        x_layer=x_layer,
        output_format="h5ad",
    )

    # Verify call
    assert len(extract_calls) == 1
    call_kwargs = extract_calls[0]
    assert call_kwargs["experiment_uri"] == experiment_uri
    assert call_kwargs["value_filter"] == value_filter
    assert call_kwargs["joinid_range"] == joinid_range
    assert call_kwargs["output_path"] == output_path
    assert call_kwargs["obs_columns"] == obs_columns
    assert call_kwargs["var_columns"] == var_columns
    assert call_kwargs["x_layer"] == x_layer

    # Verify return value
    assert result == (idx, str(output_path))


def test_extract_ranges_happy_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify extract_ranges processes all ranges in parallel.
    """
    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter='tissue == "lung"',
        id_ranges=[
            IdContiguousRange(start=0, end=10),
            IdContiguousRange(start=11, end=20),
        ],
        total_cells=20,
        range_size=10,
        output_chunk_size=10,
        num_output_chunks=2,
        output_chunk_indexes=[0, 1],
        filters={"tissue__eq": "lung"},
        obs_columns=["cell_type"],
        var_columns=["symbol"],
        x_layer="raw",
    )

    output_dir = tmp_path / "output"

    # Mock _extract_range_worker
    def fake_worker(
        idx: int,
        experiment_uri: str,
        value_filter: str,
        joinid_range: IdContiguousRange,
        output_path: Path,
        obs_columns: list[str] | None,
        var_columns: list[str] | None,
        x_layer: str,
        output_format: str,
    ) -> tuple[int, str]:
        # Create the output file/directory
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_format == "zarr":
            output_path.mkdir(exist_ok=True)
        else:
            output_path.touch()
        return idx, str(output_path)

    monkeypatch.setattr(extract, "_extract_range_worker", fake_worker)

    # Mock ProcessPoolExecutor to run synchronously
    executor = FakeExecutor(max_workers=2, worker_fn=fake_worker)
    monkeypatch.setattr(extract, "ProcessPoolExecutor", lambda **kwargs: executor)

    # Mock as_completed to return futures
    def fake_as_completed(futures: dict[object, int]) -> list[object]:
        return list(futures.keys())

    monkeypatch.setattr(extract, "as_completed", fake_as_completed)

    # Execute
    extract.extract_ranges(
        curriculum_metadata=curriculum_metadata,
        output_dir=output_dir,
        max_workers=2,
    )

    # Verify output directory was created
    assert output_dir.exists()

    # Verify worker was called for each range
    assert len(executor.calls) == 2

    # Verify output files were created (default format is zarr)
    assert (output_dir / "range_000000.zarr").exists()
    assert (output_dir / "range_000001.zarr").exists()


def test_extract_ranges_failure_handling(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify extract_ranges raises SomaExtractError when workers fail.
    """
    curriculum_metadata = SomaCurriculumMetadata(
        experiment_uri="gs://bucket/soma",
        value_filter="",
        id_ranges=[
            IdContiguousRange(start=0, end=10),
            IdContiguousRange(start=11, end=20),
        ],
        total_cells=20,
        range_size=10,
        output_chunk_size=10,
        num_output_chunks=2,
        output_chunk_indexes=[0, 1],
        filters=None,
    )

    output_dir = tmp_path / "output"

    # Mock ProcessPoolExecutor with failing worker
    class FailingExecutor(FakeExecutor):
        def submit(self, fn: object, *args: object, **kwargs: object) -> FakeFuture:
            idx = args[0]  # type: ignore[index]
            # Make second worker fail
            if idx == 1:
                future = FakeFuture(exception=RuntimeError(f"Worker {idx} failed"))
            else:
                future = FakeFuture(result=(idx, f"path_{idx}"))
            self.futures.append(future)
            return future

    executor = FailingExecutor(max_workers=2)
    monkeypatch.setattr(extract, "ProcessPoolExecutor", lambda **kwargs: executor)

    # Mock as_completed to return futures
    def fake_as_completed(futures: dict[object, int]) -> list[object]:
        return list(futures.keys())

    monkeypatch.setattr(extract, "as_completed", fake_as_completed)

    # Execute and expect error
    with pytest.raises(SomaExtractError):
        extract.extract_ranges(
            curriculum_metadata=curriculum_metadata,
            output_dir=output_dir,
        )


def test_shuffle_extracted_chunks_invalid_chunk_size(tmp_path: Path) -> None:
    """
    Verify invalid chunk_size raises ValueError.
    """
    with pytest.raises(ValueError):
        extract.shuffle_extracted_chunks(
            input_dir=tmp_path,
            output_dir=tmp_path,
            chunk_size=0,
        )


def test_shuffle_extracted_chunks_invalid_output_format(tmp_path: Path) -> None:
    """
    Verify invalid output_format raises ValueError.
    """
    with pytest.raises(ValueError):
        extract.shuffle_extracted_chunks(
            input_dir=tmp_path,
            output_dir=tmp_path,
            chunk_size=10,
            output_format="invalid",
        )


def test_shuffle_extracted_chunks_no_input_files(tmp_path: Path) -> None:
    """
    Verify missing input files raises ValueError.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    with pytest.raises(ValueError):
        extract.shuffle_extracted_chunks(
            input_dir=input_dir,
            output_dir=output_dir,
            chunk_size=10,
        )


def test_shuffle_extracted_chunks_happy_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """
    Verify shuffle_extracted_chunks redistributes cells correctly.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"

    # Create fake input files
    (input_dir / "range_000000.h5ad").touch()
    (input_dir / "range_000001.h5ad").touch()

    # Mock anndata.read_h5ad
    def fake_read_h5ad(path: Path, backed: str | None = None) -> FakeAnnData:
        if backed:
            assert backed == "r"
        # Return fake with 3 cells per file
        return FakeAnnData(n_obs=3)

    monkeypatch.setattr("anndata.read_h5ad", fake_read_h5ad)

    # Mock np.random.permutation
    def fake_permutation(n: int) -> np.ndarray:
        # Return deterministic permutation
        return np.arange(n)

    monkeypatch.setattr(np.random, "permutation", fake_permutation)

    # Mock _write_shuffled_chunk worker
    def fake_write_chunk(
        chunk_idx: int,
        output_chunk_idx: int,
        chunk_indices: object,
        input_files: object,
        input_format: str,
        output_dir: Path,
        output_format: str,
        var_joinids: list[int] | None,
    ) -> tuple[int, int, str]:
        # Create output file
        output_dir.mkdir(parents=True, exist_ok=True)
        if output_format == "zarr":
            output_path = output_dir / f"extract_{output_chunk_idx:06d}.zarr"
            output_path.mkdir(exist_ok=True)
        else:
            output_path = output_dir / f"extract_{output_chunk_idx:06d}.h5ad"
            output_path.touch()
        return chunk_idx, output_chunk_idx, str(output_path)

    # Mock ProcessPoolExecutor
    executor = FakeExecutor(max_workers=1, worker_fn=fake_write_chunk)
    monkeypatch.setattr(extract, "ProcessPoolExecutor", lambda **kwargs: executor)

    # Mock as_completed
    def fake_as_completed(futures: dict[object, int]) -> list[object]:
        return list(futures.keys())

    monkeypatch.setattr(extract, "as_completed", fake_as_completed)

    # Execute
    extract.shuffle_extracted_chunks(
        input_dir=input_dir,
        output_dir=output_dir,
        chunk_size=4,
        input_format="h5ad",
        output_format="h5ad",
        max_workers=1,
    )

    # Verify output directory was created
    assert output_dir.exists()

    # With 6 total cells and chunk_size=4, should create 2 extracts
    assert (output_dir / "extract_000000.h5ad").exists()
    assert (output_dir / "extract_000001.h5ad").exists()
