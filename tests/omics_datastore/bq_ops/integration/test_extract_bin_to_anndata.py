from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from cellarium.nexus.omics_datastore.bq_ops.extract import extract as extract_module
from cellarium.nexus.shared import schemas as schemas_module


def test_extract_bin_to_anndata_end_to_end(tmp_path: Path) -> None:
    """
    Exercise DataExtractor.extract_bin_to_anndata end-to-end using stubbed BigQuery calls.

    :param tmp_path: Temporary directory to write the output .h5ad file
    """
    # Build a DataExtractor with dummy client; monkeypatch method surfaces
    de = extract_module.DataExtractor(client=object(), project="P", dataset="D", extract_table_prefix="X_")

    # Fake implementations for BigQuery calls
    de.get_features = lambda: [
        schemas_module.FeatureSchema(id=1, symbol="A", ensemble_id="ensA"),
        schemas_module.FeatureSchema(id=2, symbol="B", ensemble_id="ensB"),
        schemas_module.FeatureSchema(id=3, symbol="C", ensemble_id="ensC"),
    ]  # type: ignore[assignment]

    de.get_cells_in_bin_range = lambda start_bin, end_bin, obs_columns=None: [
        {"id": 1, "ingest_id": "ing1", "group": "g1", "score": 1.0},
        {"id": 2, "ingest_id": "ing2", "group": None, "score": None},
    ]  # type: ignore[assignment]

    de.get_matrix_data = lambda start_bin, end_bin: [
        {"cell_id": 1, "feature_data": [{"feature_id": 1, "raw_counts": 1.0}, {"feature_id": 3, "raw_counts": 2.0}]},
        {"cell_id": 2, "feature_data": [{"feature_id": 2, "raw_counts": 3.5}]},
    ]  # type: ignore[assignment]

    meta = schemas_module.ExtractMetadata(
        total_bins=1,
        last_bin_size=2,
        category_metadata={"group": ["g1", "g2"]},
    )

    out = tmp_path / "b5.h5ad"
    de.extract_bin_to_anndata(
        bin_number=5,
        output_path=out,
        extract_metadata=meta,
        obs_columns=["group", "score"],
    )

    assert out.exists()
    adata = ad.read_h5ad(filename=out)

    assert adata.shape == (2, 3)
    assert list(adata.obs.index) == ["1", "2"]
    assert list(adata.var.index) == ["ensA", "ensB", "ensC"]
    assert isinstance(adata.obs["group"].dtype, pd.CategoricalDtype)  # type: ignore[arg-type]
    assert list(adata.obs["group"].cat.categories) == ["g1", "g2"]  # type: ignore[attr-defined]
    # missing string in categorical column becomes NaN after categorical conversion
    assert pd.isna(adata.obs.loc["2", "group"])  # type: ignore[arg-type]
    # missing numeric -> NaN
    assert np.isnan(adata.obs.loc["2", "score"])  # type: ignore[arg-type]
    # nnz count
    assert adata.X.nnz == 3
