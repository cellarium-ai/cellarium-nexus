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

    # Create schema matching the adata
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string")],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=pd.DataFrame(index=adata.var.index)),
        x_validation_type="count_matrix",
    )

    sanitize_for_ingest(adata=adata, ingest_schema=schema)

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
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string", nullable=False)],
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

    # Schema defines canonical var - no columns, only feature IDs
    schema_var_df = pd.DataFrame(index=all_features)
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=schema_var_df),
        x_validation_type="count_matrix",
    )

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
    sanitize_for_ingest(adata=adata1, ingest_schema=schema)
    sanitize_for_ingest(adata=adata2, ingest_schema=schema)
    h5ad_path1 = tmp_path / "file1.h5ad"
    h5ad_path2 = tmp_path / "file2.h5ad"
    adata1.write_h5ad(filename=h5ad_path1)
    adata2.write_h5ad(filename=h5ad_path2)

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


def test_sanitize_all_supported_types_write_h5ad(tmp_path) -> None:
    """
    Verify comprehensive type support for sanitization and H5AD writing.

    Tests all supported obs column types defined in SomaDtype in both
    nullable and non-nullable configurations (where applicable), verifying they can be:
    1. Sanitized correctly (mapped to correct dtype)
    2. Written to .h5ad without error
    """
    n_obs = 5
    adata = anndata.AnnData(X=sp.csr_matrix((n_obs, 2)))

    # Helper to create range data
    def make_range(dtype, n=n_obs):
        return np.arange(n).astype(dtype)

    # Generate data covering all types
    # Note: Use lists for nullable columns to avoid index alignment issues during temporary DataFrame construction
    # (pd.Series would require explicit index matching "c0".."c4")
    obs_data = {
        # Boolean
        "bool_req": [True, False, True, False, True],
        "bool_null": [True, False, None, True, False],
        # Signed integers
        "int8_req": make_range("int8"),
        "int8_null": [1, 2, None, 4, 5],
        "int16_req": make_range("int16"),
        "int16_null": [1, 2, None, 4, 5],
        "int32_req": make_range("int32"),
        "int32_null": [1, 2, None, 4, 5],
        "int64_req": make_range("int64"),
        "int64_null": [1, 2, None, 4, 5],
        # Floating point
        "float32_req": np.array([1.1, 2.2, 3.3, 4.4, 5.5], dtype="float32"),
        "float32_null": [1.1, 2.2, None, 4.4, np.nan],
        "float64_req": np.array([1.1, 2.2, 3.3, 4.4, 5.5], dtype="float64"),
        "float64_null": [1.1, 2.2, None, 4.4, np.nan],
        # String types
        "str_req": ["a", "b", "c", "d", "e"],
        "str_null": ["a", "b", None, "d", None],
        # Categorical (SOMA treats as Enum, here we test basic category IO)
        "cat_req": pd.Categorical(["a", "b", "a", "b", "c"]),
        # Nullable categoricals are supported by pandas but SOMA schema maps to 'category'
        "cat_null": pd.Categorical(["a", None, "b", None, "c"]),
    }

    adata.obs = pd.DataFrame(obs_data, index=[f"c{i}" for i in range(n_obs)])

    obs_columns = [
        # Bool
        ObsSchemaDescriptor(name="bool_req", dtype="bool", nullable=False),
        ObsSchemaDescriptor(name="bool_null", dtype="bool", nullable=True),
        # Ints
        ObsSchemaDescriptor(name="int8_req", dtype="int8", nullable=False),
        ObsSchemaDescriptor(name="int8_null", dtype="int8", nullable=True),
        ObsSchemaDescriptor(name="int16_req", dtype="int16", nullable=False),
        ObsSchemaDescriptor(name="int16_null", dtype="int16", nullable=True),
        ObsSchemaDescriptor(name="int32_req", dtype="int32", nullable=False),
        ObsSchemaDescriptor(name="int32_null", dtype="int32", nullable=True),
        ObsSchemaDescriptor(name="int64_req", dtype="int64", nullable=False),
        ObsSchemaDescriptor(name="int64_null", dtype="int64", nullable=True),
        # Floats
        ObsSchemaDescriptor(name="float32_req", dtype="float32", nullable=False),
        ObsSchemaDescriptor(name="float32_null", dtype="float32", nullable=True),
        ObsSchemaDescriptor(name="float64_req", dtype="float64", nullable=False),
        ObsSchemaDescriptor(name="float64_null", dtype="float64", nullable=True),
        # String
        ObsSchemaDescriptor(name="str_req", dtype="string", nullable=False),
        ObsSchemaDescriptor(name="str_null", dtype="string", nullable=True),
        # Category
        ObsSchemaDescriptor(name="cat_req", dtype="category", nullable=False),
        ObsSchemaDescriptor(name="cat_null", dtype="category", nullable=True),
    ]

    schema = IngestSchema(
        obs_columns=obs_columns,
        var_schema=ExperimentVarSchema.from_dataframe(var_df=pd.DataFrame(index=adata.var.index)),
        x_validation_type="count_matrix",
    )

    sanitize_for_ingest(adata=adata, ingest_schema=schema)

    # Verify Ints (Standard Nullable Mapping)
    assert str(adata.obs["int32_null"].dtype) == "Int32"
    assert adata.obs["int32_null"].isna().sum() == 1

    # Verify Bools (Standard Nullable Mapping)
    assert str(adata.obs["bool_null"].dtype) == "boolean"
    assert adata.obs["bool_null"].isna().sum() == 1

    # Verify Floats (Special Mapping: float32/64 to support H5AD write)
    assert adata.obs["float32_null"].dtype == "float32"
    assert adata.obs["float32_null"].isna().sum() == 2  # None and np.nan both become nan
    assert adata.obs["float64_null"].dtype == "float64"
    assert adata.obs["float64_null"].isna().sum() == 2

    # Verify Strings (Special Mapping: object to support H5AD write)
    assert adata.obs["str_null"].dtype == "object"
    assert adata.obs["str_null"].isna().sum() == 0  # None values are sanitized to "unknown"
    assert (adata.obs["str_null"] == "unknown").sum() == 2  # Two None values became "unknown"

    # Verify Categories
    assert isinstance(adata.obs["cat_req"].dtype, pd.CategoricalDtype)
    assert isinstance(adata.obs["cat_null"].dtype, pd.CategoricalDtype)
    assert adata.obs["cat_null"].isna().sum() == 2

    # Verify Write
    out_path = tmp_path / "test_full_types.h5ad"

    # This should not raise IORegistryError
    adata.write_h5ad(out_path)

    # Verify Read Back matches
    adata_read = anndata.read_h5ad(out_path)

    # Content integrity check samples
    assert adata_read.obs["str_req"].tolist() == ["a", "b", "c", "d", "e"]
    assert adata_read.obs["str_null"].iloc[0] == "a"
    assert adata_read.obs["int32_req"].iloc[0] == 0
    assert adata_read.obs["float64_req"].iloc[0] == 1.1
