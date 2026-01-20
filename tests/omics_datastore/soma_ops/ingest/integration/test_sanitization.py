"""
Integration tests for SOMA ingest sanitization.

Test sanitize_for_ingest and sanitize_first_adata_for_schema
to verify sanitized AnnData is compatible with SOMA registration.
"""

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledbsoma.io

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import (
    sanitize_first_adata_for_schema,
    sanitize_for_ingest,
)
from cellarium.nexus.shared.schemas.omics_datastore import (
    ExperimentVarSchema,
    IngestSchema,
    ObsSchemaDescriptor,
)

# Tests for sanitize_for_ingest


def test_sanitize_for_ingest_produces_soma_compatible_adata(
    small_valid_anndata: tuple[anndata.AnnData, str],
    soma_experiment_uri: str,
) -> None:
    """Verify sanitized AnnData can be ingested into SOMA experiment."""
    adata, _ = small_valid_anndata

    # Add extra slots that should be removed
    adata.uns["metadata"] = {"key": "value"}
    adata.obsm["X_pca"] = np.random.rand(adata.n_obs, 10)
    adata.layers["counts"] = adata.X.copy()  # type: ignore[union-attr]

    sanitize_for_ingest(adata=adata)

    # Should be able to create SOMA experiment from sanitized adata
    tiledbsoma.io.from_anndata(
        experiment_uri=soma_experiment_uri,
        anndata=adata,
        measurement_name="RNA",
        obs_id_name="obs_id",
        var_id_name="var_id",
        ingest_mode="schema_only",
    )

    # Verify experiment was created
    with tiledbsoma.Experiment.open(uri=soma_experiment_uri, mode="r") as exp:
        assert "RNA" in exp.ms
        assert "obs" in exp
        assert "var" in exp.ms["RNA"]


# Tests for sanitize_first_adata_for_schema


def test_sanitize_first_adata_for_schema_replaces_var_from_schema(
    tmp_path,
) -> None:
    """Verify var is replaced entirely from schema var index (no columns)."""
    n_obs = 5

    X = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0], [0, 0, 5], [1, 1, 0]], dtype=np.int32))

    obs = pd.DataFrame(
        data={"cell_type": ["a", "b", "c", "d", "e"]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        data={"feature_name": ["gene_0", "gene_1", "gene_2"]},
        index=["ENSG0001", "ENSG0002", "ENSG0003"],
    )

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Schema has more features than adata - var will be replaced with schema index only
    schema_var_df = pd.DataFrame(
        index=["ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004", "ENSG0005"],
    )
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=schema_var_df),
        x_validation_type="count_matrix",
    )

    sanitize_first_adata_for_schema(adata=adata, ingest_schema=schema)

    # Verify var now has all 5 features from schema (no columns)
    assert adata.n_vars == 5
    assert list(adata.var.index) == ["ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004", "ENSG0005"]
    assert len(adata.var.columns) == 0

    # Verify X shape expanded
    assert adata.X.shape == (5, 5)  # type: ignore[union-attr]

    # Verify original data preserved in first 3 columns
    X_dense = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X  # type: ignore[union-attr]
    assert X_dense[0, 0] == 1  # type: ignore[index]
    assert X_dense[0, 2] == 2  # type: ignore[index]

    # Verify new columns are zero-filled
    assert X_dense[:, 3].sum() == 0  # type: ignore[index]
    assert X_dense[:, 4].sum() == 0  # type: ignore[index]


def test_sanitize_first_adata_for_schema_enables_multi_file_registration(
    tmp_path,
) -> None:
    """Verify sanitized first adata can be used for multi-file registration."""
    n_obs = 5
    all_features = ["ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004", "ENSG0005"]

    # Create first adata with subset of features
    X1 = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0], [0, 0, 5], [1, 1, 0]], dtype=np.int32))
    obs1 = pd.DataFrame(data={"cell_type": ["a", "b", "c", "d", "e"]}, index=[f"cell1_{i}" for i in range(n_obs)])
    var1 = pd.DataFrame(index=["ENSG0001", "ENSG0002", "ENSG0003"])  # No columns
    adata1 = anndata.AnnData(X=X1, obs=obs1, var=var1)

    # Create second adata with different subset of features
    X2 = sp.csr_matrix(np.array([[2, 1], [0, 4], [3, 0]], dtype=np.int32))
    obs2 = pd.DataFrame(data={"cell_type": ["x", "y", "z"]}, index=[f"cell2_{i}" for i in range(3)])
    var2 = pd.DataFrame(index=["ENSG0004", "ENSG0005"])  # No columns
    adata2 = anndata.AnnData(X=X2, obs=obs2, var=var2)

    # Sanitize and save to h5ad
    sanitize_for_ingest(adata=adata1)
    sanitize_for_ingest(adata=adata2)
    h5ad_path1 = tmp_path / "file1.h5ad"
    h5ad_path2 = tmp_path / "file2.h5ad"
    adata1.write_h5ad(filename=h5ad_path1)
    adata2.write_h5ad(filename=h5ad_path2)

    # Schema defines canonical var - no columns, only feature IDs
    schema_var_df = pd.DataFrame(index=all_features)
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=schema_var_df),
        x_validation_type="count_matrix",
    )

    # Load and sanitize first adata for schema creation
    first_adata = anndata.read_h5ad(filename=h5ad_path1)
    sanitize_first_adata_for_schema(adata=first_adata, ingest_schema=schema)

    # Verify first adata has all features now (from schema)
    assert first_adata.n_vars == 5
    assert list(first_adata.var.index) == all_features
    assert len(first_adata.var.columns) == 0

    # Create SOMA experiment from sanitized first adata
    experiment_uri = str(tmp_path / "experiment")
    tiledbsoma.io.from_anndata(
        experiment_uri=experiment_uri,
        anndata=adata1,
        measurement_name="RNA",
        obs_id_name="obs_id",
        var_id_name="var_id",
        ingest_mode="schema_only",
    )

    # Register both h5ad files
    registration_mapping = tiledbsoma.io.register_h5ads(
        experiment_uri=experiment_uri,
        h5ad_file_names=[str(h5ad_path1), str(h5ad_path2)],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    # Verify registration succeeded
    assert registration_mapping is not None
