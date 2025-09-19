from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import extract as extract_module
from cellarium.nexus.shared import schemas as schemas_module


def test_convert_matrix_data_to_coo_happy_path() -> None:
    """
    Validate that convert_matrix_data_to_coo produces expected row/col/data lists.
    """
    matrix_data = [
        {"cell_id": 10, "feature_data": [{"feature_id": 1, "raw_counts": 2.0}, {"feature_id": 2, "raw_counts": 0.0}]},
        {"cell_id": 11, "feature_data": [{"feature_id": 2, "raw_counts": 3.5}]},
    ]
    cell_index_to_row_num = {10: 0, 11: 1}
    feature_id_to_col_num = {1: 0, 2: 1}

    rows, cols, data = extract_module.DataExtractor.convert_matrix_data_to_coo(
        self=None,  # type: ignore[arg-type]
        matrix_data=matrix_data,
        cell_index_to_row_num=cell_index_to_row_num,
        feature_id_to_col_num=feature_id_to_col_num,
    )

    assert rows == [0, 0, 1]
    assert cols == [0, 1, 1]
    assert data == [2.0, 0.0, 3.5]


def test_convert_matrix_data_to_coo_raises_on_invalid_ids() -> None:
    """
    Ensure that convert_matrix_data_to_coo raises when cell or feature IDs are unmapped.
    """
    # unknown cell id
    matrix_data = [{"cell_id": 999, "feature_data": [{"feature_id": 1, "raw_counts": 1.0}]}]
    with pytest.raises(KeyError):
        extract_module.DataExtractor.convert_matrix_data_to_coo(
            self=None,  # type: ignore[arg-type]
            matrix_data=matrix_data,
            cell_index_to_row_num={10: 0},
            feature_id_to_col_num={1: 0},
        )

    # unknown feature id
    matrix_data2 = [{"cell_id": 10, "feature_data": [{"feature_id": 999, "raw_counts": 1.0}]}]
    with pytest.raises(KeyError):
        extract_module.DataExtractor.convert_matrix_data_to_coo(
            self=None,  # type: ignore[arg-type]
            matrix_data=matrix_data2,
            cell_index_to_row_num={10: 0},
            feature_id_to_col_num={1: 0},
        )


def test_extract_bin_to_anndata_writes_h5ad_with_expected_obs_var(tmp_path: Path) -> None:
    """
    Exercise extract_bin_to_anndata to write a valid .h5ad with correct obs/var and metadata handling.

    :param tmp_path: Temporary directory for output file
    """
    # Build extractor with dummy client
    de = extract_module.DataExtractor(client=object(), project="p", dataset="d", extract_table_prefix="x_")

    # Features: ids map to var (ensemble ids)
    features = [
        schemas_module.FeatureSchema(id=1, symbol="A", ensemble_id="ensA"),
        schemas_module.FeatureSchema(id=2, symbol="B", ensemble_id="ensB"),
    ]

    cells = [
        {"id": 100, "ingest_id": "ing1", "sex": "M", "score": 1.5, "label": None, "all_missing": None},
        {"id": 101, "ingest_id": "ing2", "sex": None, "score": None, "label": "x", "all_missing": None},
    ]

    matrix_rows = [
        {"cell_id": 100, "feature_data": [{"feature_id": 1, "raw_counts": 2.0}]},
        {"cell_id": 101, "feature_data": [{"feature_id": 2, "raw_counts": 3.0}]},
    ]

    # Monkeypatch extractor methods
    de.get_features = lambda: features  # type: ignore[assignment]
    de.get_cells_in_bin_range = lambda start_bin, end_bin, obs_columns=None: cells  # type: ignore[assignment]
    de.get_matrix_data = lambda start_bin, end_bin: matrix_rows  # type: ignore[assignment]

    meta = schemas_module.ExtractMetadata(
        total_bins=1,
        last_bin_size=2,
        category_metadata={"sex": ["M", "F"]},
    )

    out_path = tmp_path / "out.h5ad"
    de.extract_bin_to_anndata(
        bin_number=5,
        output_path=out_path,
        extract_metadata=meta,
        obs_columns=["sex", "score", "label", "all_missing"],
    )

    assert out_path.exists()

    adata = ad.read_h5ad(filename=out_path)
    # shapes
    assert adata.shape == (2, 2)
    # var index is ensemble ids
    assert list(adata.var.index) == ["ensA", "ensB"]
    # obs index is original ids (coerced to str by AnnData I/O)
    assert list(adata.obs.index) == ["100", "101"]

    # categorical
    assert isinstance(adata.obs["sex"].dtype, pd.CategoricalDtype)  # type: ignore[arg-type]
    assert list(adata.obs["sex"].cat.categories) == ["M", "F"]  # type: ignore[attr-defined]

    # string None becomes empty
    assert list(adata.obs["label"]) == ["", "x"]

    # numeric None becomes NaN
    assert np.isnan(adata.obs.loc["101", "score"])  # type: ignore[arg-type]

    # all-None column skipped
    assert "all_missing" not in adata.obs.columns

    # dtype and nnz
    assert adata.X.dtype == np.float32
    assert adata.X.nnz == 2
