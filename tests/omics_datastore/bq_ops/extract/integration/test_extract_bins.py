import typing
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import extract as extract_module
from cellarium.nexus.shared import schemas as schemas_module


def _write_min_h5ad(path: Path) -> None:
    X = np.zeros((1, 1), dtype=np.float32)
    obs = pd.DataFrame(index=["c1"])  # index dtype string ok
    var = pd.DataFrame(index=["g1"])  # feature id
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(filename=path)


def test_extract_bins_happy_path(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, patched_process_pool_executor: typing.Iterator[None]
) -> None:
    """
    Exercise extract_bins happy path by patching worker to write minimal .h5ad outputs.

    :param monkeypatch: Pytest monkeypatch fixture
    :param tmp_path: Temporary directory to write output files
    :param patched_process_pool_executor: Fixture patching ProcessPoolExecutor to threads
    """

    # Patch worker to write small files synchronously
    def _worker(**kwargs: typing.Any) -> None:
        out: Path = kwargs["output_path"]
        _write_min_h5ad(out)

    monkeypatch.setattr(extract_module, "extract_bin_to_anndata_worker", _worker)

    meta = schemas_module.ExtractMetadata(total_bins=2, last_bin_size=1)

    extract_module.extract_bins(
        client=object(),
        project="P",
        dataset="D",
        extract_table_prefix="X_",
        bins=[10, 11],
        output_dir=tmp_path,
        extract_metadata=meta,
        obs_columns=None,
        max_workers=2,
    )

    p10 = tmp_path / "extract_10.h5ad"
    p11 = tmp_path / "extract_11.h5ad"
    assert p10.exists() and p11.exists()
    # sanity read
    ad.read_h5ad(filename=p10)
    ad.read_h5ad(filename=p11)


def test_extract_bins_surfaces_worker_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, patched_process_pool_executor: typing.Iterator[None]
) -> None:
    """
    Verify that extract_bins propagates exceptions raised by extract_bin_to_anndata_worker.

    :param monkeypatch: Pytest monkeypatch fixture
    :param tmp_path: Temporary directory for outputs
    :param patched_process_pool_executor: Fixture patching ProcessPoolExecutor to threads
    """

    def _worker(**kwargs: typing.Any) -> None:
        raise ValueError("fail")

    monkeypatch.setattr(extract_module, "extract_bin_to_anndata_worker", _worker)

    meta = schemas_module.ExtractMetadata(total_bins=1, last_bin_size=1)

    with pytest.raises(ValueError):
        extract_module.extract_bins(
            client=object(),
            project="P",
            dataset="D",
            extract_table_prefix="X_",
            bins=[5],
            output_dir=tmp_path,
            extract_metadata=meta,
            obs_columns=None,
            max_workers=1,
        )
