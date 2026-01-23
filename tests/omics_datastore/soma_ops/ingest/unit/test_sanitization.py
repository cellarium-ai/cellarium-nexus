import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import sanitization
from cellarium.nexus.shared.schemas.omics_datastore import (
    ExperimentVarSchema,
    IngestSchema,
    ObsSchemaDescriptor,
)


@pytest.fixture()
def adata_with_all_slots() -> anndata.AnnData:
    """
    Create an AnnData with all slots populated for sanitization testing.

    :return: AnnData with obsm, varm, uns, obsp, varp, layers, and raw.
    """
    n_obs, n_vars = 5, 4
    X = sp.random(n_obs, n_vars, density=0.3, format="csr", dtype=np.float32)

    obs = pd.DataFrame(
        data={"cell_type": ["A", "B", "C", "D", "E"]},
        index=pd.Index([f"cell_{i}" for i in range(n_obs)], name="obs_index"),
    )
    var = pd.DataFrame(
        data={"gene_name": [f"gene_{j}" for j in range(n_vars)]},
        index=pd.Index([f"ENSG{j:04d}" for j in range(n_vars)], name="var_index"),
    )

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Populate obsm
    adata.obsm["X_pca"] = np.random.rand(n_obs, 10).astype(np.float32)
    adata.obsm["X_umap"] = np.random.rand(n_obs, 2).astype(np.float32)

    # Populate varm
    adata.varm["PCs"] = np.random.rand(n_vars, 10).astype(np.float32)

    # Populate uns
    adata.uns["title"] = "test_dataset"
    adata.uns["neighbors"] = {"connectivities_key": "connectivities"}

    # Populate obsp
    adata.obsp["distances"] = sp.random(n_obs, n_obs, density=0.2, format="csr")
    adata.obsp["connectivities"] = sp.random(n_obs, n_obs, density=0.2, format="csr")

    # Populate varp
    adata.varp["gene_correlations"] = sp.random(n_vars, n_vars, density=0.3, format="csr")

    # Populate layers
    adata.layers["counts"] = sp.random(n_obs, n_vars, density=0.3, format="csr", dtype=np.int32)
    adata.layers["normalized"] = sp.random(n_obs, n_vars, density=0.3, format="csr", dtype=np.float32)

    # Set raw
    adata.raw = adata.copy()

    return adata


# Tests for _remove_obsm


def test_remove_obsm_removes_all_keys(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all obsm keys are removed."""
    assert len(adata_with_all_slots.obsm) > 0

    sanitization._remove_obsm(adata=adata_with_all_slots)

    assert len(adata_with_all_slots.obsm) == 0


def test_remove_obsm_handles_empty() -> None:
    """Verify function handles AnnData with no obsm."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization._remove_obsm(adata=adata)

    assert len(adata.obsm) == 0


# Tests for _remove_varm


def test_remove_varm_removes_all_keys(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all varm keys are removed."""
    assert len(adata_with_all_slots.varm) > 0

    sanitization._remove_varm(adata=adata_with_all_slots)

    assert len(adata_with_all_slots.varm) == 0


def test_remove_varm_handles_empty() -> None:
    """Verify function handles AnnData with no varm."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization._remove_varm(adata=adata)

    assert len(adata.varm) == 0


# Tests for _remove_uns


def test_remove_uns_removes_all_keys(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all uns keys are removed."""
    assert len(adata_with_all_slots.uns) > 0

    sanitization._remove_uns(adata=adata_with_all_slots)

    assert len(adata_with_all_slots.uns) == 0


def test_remove_uns_handles_empty() -> None:
    """Verify function handles AnnData with no uns."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization._remove_uns(adata=adata)

    assert len(adata.uns) == 0


# Tests for _remove_obsp


def test_remove_obsp_removes_all_keys(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all obsp keys are removed."""
    assert len(adata_with_all_slots.obsp) > 0

    sanitization._remove_obsp(adata=adata_with_all_slots)

    assert len(adata_with_all_slots.obsp) == 0


def test_remove_obsp_handles_empty() -> None:
    """Verify function handles AnnData with no obsp."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization._remove_obsp(adata=adata)

    assert len(adata.obsp) == 0


# Tests for _remove_varp


def test_remove_varp_removes_all_keys(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all varp keys are removed."""
    assert len(adata_with_all_slots.varp) > 0

    sanitization._remove_varp(adata=adata_with_all_slots)

    assert len(adata_with_all_slots.varp) == 0


def test_remove_varp_handles_empty() -> None:
    """Verify function handles AnnData with no varp."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization._remove_varp(adata=adata)

    assert len(adata.varp) == 0


# Tests for _remove_layers


def test_remove_layers_removes_all(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all layers are removed."""
    assert len(adata_with_all_slots.layers) > 0

    sanitization._remove_layers(adata=adata_with_all_slots)

    assert len(adata_with_all_slots.layers) == 0


def test_remove_layers_handles_empty() -> None:
    """Verify function handles AnnData with no layers."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization._remove_layers(adata=adata)

    assert len(adata.layers) == 0


# Tests for _remove_raw


def test_remove_raw_removes(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify raw is removed."""
    assert adata_with_all_slots.raw is not None

    sanitization._remove_raw(adata=adata_with_all_slots)

    assert adata_with_all_slots.raw is None


def test_remove_raw_handles_no_raw() -> None:
    """Verify function handles AnnData with no raw."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    assert adata.raw is None

    sanitization._remove_raw(adata=adata)

    assert adata.raw is None


# Tests for _reset_index_names


def test_reset_index_names_resets_both(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify obs and var index names are reset to None."""
    assert adata_with_all_slots.obs.index.name == "obs_index"
    assert adata_with_all_slots.var.index.name == "var_index"

    sanitization._reset_index_names(adata=adata_with_all_slots)

    assert adata_with_all_slots.obs.index.name is None
    assert adata_with_all_slots.var.index.name is None


def test_reset_index_names_handles_no_names() -> None:
    """Verify function handles AnnData with no index names."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    assert adata.obs.index.name is None
    assert adata.var.index.name is None

    sanitization._reset_index_names(adata=adata)

    assert adata.obs.index.name is None
    assert adata.var.index.name is None


@pytest.fixture()
def basic_ingest_schema():
    """Create an IngestSchema matching adata_with_all_slots fixture."""

    # Create schema var DataFrame matching adata_with_all_slots (4 features)
    var_df = pd.DataFrame(index=[f"ENSG{j:04d}" for j in range(4)])
    return IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string")],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )


# Tests for _sanitize_var_metadata_with_schema


def test_sanitize_var_metadata_with_schema_updates_columns() -> None:
    """Verify var columns are updated from schema, index preserved."""

    # Create schema with extra columns
    var_df = pd.DataFrame(data={"gene_symbol": ["A", "B", "C", "D"]}, index=[f"ENSG{j:04d}" for j in range(4)])
    ingest_schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string")],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    # Adata with different columns
    adata = anndata.AnnData(
        X=sp.csr_matrix((5, 4)),
        var=pd.DataFrame(data={"old_col": range(4)}, index=[f"ENSG{j:04d}" for j in range(4)]),
    )

    sanitization._sanitize_var_metadata_with_schema(adata=adata, ingest_schema=ingest_schema)

    assert list(adata.var.columns) == ["gene_symbol"]
    assert list(adata.var["gene_symbol"]) == ["A", "B", "C", "D"]


def test_sanitize_var_metadata_with_schema_handles_subset() -> None:
    """Verify function handles AnnData with subset of features."""
    from cellarium.nexus.shared.schemas.omics_datastore import (
        ExperimentVarSchema,
        IngestSchema,
    )

    var_df = pd.DataFrame(data={"gene_symbol": ["A", "B", "C"]}, index=["g1", "g2", "g3"])
    ingest_schema = IngestSchema(
        obs_columns=[],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    adata = anndata.AnnData(
        X=sp.csr_matrix((1, 2)),
        var=pd.DataFrame(index=["g1", "g3"]),
    )

    sanitization._sanitize_var_metadata_with_schema(adata=adata, ingest_schema=ingest_schema)

    assert list(adata.var.index) == ["g1", "g3"]
    assert list(adata.var["gene_symbol"]) == ["A", "C"]


# Tests for sanitize_for_ingest


def test_sanitize_for_ingest_removes_all_unsupported_slots(
    adata_with_all_slots: anndata.AnnData, basic_ingest_schema
) -> None:
    """Verify all unsupported slots are removed and var columns stripped."""
    # Verify slots are populated before sanitization
    assert len(adata_with_all_slots.obsm) > 0
    assert len(adata_with_all_slots.varm) > 0
    assert len(adata_with_all_slots.uns) > 0
    assert len(adata_with_all_slots.obsp) > 0
    assert len(adata_with_all_slots.varp) > 0
    assert len(adata_with_all_slots.layers) > 0
    assert adata_with_all_slots.raw is not None
    assert adata_with_all_slots.obs.index.name is not None
    assert adata_with_all_slots.var.index.name is not None
    assert len(adata_with_all_slots.var.columns) > 0

    sanitization.sanitize_for_ingest(adata=adata_with_all_slots, ingest_schema=basic_ingest_schema)

    # Verify all unsupported slots are removed
    assert len(adata_with_all_slots.obsm) == 0
    assert len(adata_with_all_slots.varm) == 0
    assert len(adata_with_all_slots.uns) == 0
    assert len(adata_with_all_slots.obsp) == 0
    assert len(adata_with_all_slots.varp) == 0
    assert len(adata_with_all_slots.layers) == 0
    assert adata_with_all_slots.raw is None
    assert adata_with_all_slots.obs.index.name is None
    assert adata_with_all_slots.var.index.name is None
    # Verify var columns are stripped
    assert len(adata_with_all_slots.var.columns) == 0


def test_sanitize_for_ingest_preserves_x_and_filters_obs(
    adata_with_all_slots: anndata.AnnData, basic_ingest_schema
) -> None:
    """Verify X is preserved and obs is filtered to schema columns."""
    original_x_shape = adata_with_all_slots.X.shape
    original_var_index = list(adata_with_all_slots.var.index)
    original_n_obs = adata_with_all_slots.n_obs
    original_n_vars = adata_with_all_slots.n_vars

    sanitization.sanitize_for_ingest(adata=adata_with_all_slots, ingest_schema=basic_ingest_schema)

    assert adata_with_all_slots.X.shape == original_x_shape
    # Obs columns should be filtered to schema columns only
    assert list(adata_with_all_slots.obs.columns) == ["cell_type"]
    # Var columns should be stripped, only index preserved
    assert list(adata_with_all_slots.var.columns) == []
    assert list(adata_with_all_slots.var.index) == original_var_index
    assert adata_with_all_slots.n_obs == original_n_obs
    assert adata_with_all_slots.n_vars == original_n_vars


def test_sanitize_for_ingest_operates_in_place(adata_with_all_slots: anndata.AnnData, basic_ingest_schema) -> None:
    """Verify sanitization operates in-place on the same object."""
    original_id = id(adata_with_all_slots)

    sanitization.sanitize_for_ingest(adata=adata_with_all_slots, ingest_schema=basic_ingest_schema)

    assert id(adata_with_all_slots) == original_id


def test_sanitize_for_ingest_handles_minimal_adata() -> None:
    """Verify sanitization handles minimal AnnData with only X."""
    from cellarium.nexus.shared.schemas.omics_datastore import (
        ExperimentVarSchema,
        IngestSchema,
    )

    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    # Create minimal schema matching the adata dimensions
    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization.sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)

    assert len(adata.obsm) == 0
    assert len(adata.varm) == 0
    assert len(adata.uns) == 0
    assert len(adata.obsp) == 0
    assert len(adata.varp) == 0
    assert len(adata.layers) == 0
    assert adata.raw is None


# Tests for _sanitize_obs_with_schema


def test_sanitize_obs_with_schema_filters_columns() -> None:
    """Verify obs columns are filtered to only those in schema."""

    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    adata.obs = pd.DataFrame(
        {"cell_type": ["A", "B", "C"], "extra_col": [1, 2, 3], "another_col": ["x", "y", "z"]},
        index=pd.Index(["c0", "c1", "c2"]),
    )

    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string")],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization._sanitize_obs_with_schema(adata=adata, ingest_schema=ingest_schema)

    assert list(adata.obs.columns) == ["cell_type"]
    assert list(adata.obs["cell_type"]) == ["A", "B", "C"]


def test_sanitize_obs_with_schema_casts_dtypes() -> None:
    """Verify obs columns are cast to target dtypes."""

    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    # Create obs with integer values that should be cast to string
    adata.obs = pd.DataFrame(
        {"count": [1, 2, 3]},
        index=pd.Index(["c0", "c1", "c2"]),
    )

    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="count", dtype="string")],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization._sanitize_obs_with_schema(adata=adata, ingest_schema=ingest_schema)

    # Accept either 'object' or Pandas 'string' dtype
    assert str(adata.obs["count"].dtype) in ("object", "string")
    assert list(adata.obs["count"]) == ["1", "2", "3"]


def test_sanitize_obs_with_schema_fills_nullable_missing_columns() -> None:
    """Verify missing nullable columns are filled with NA."""

    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    adata.obs = pd.DataFrame(
        {"cell_type": ["A", "B", "C"]},
        index=pd.Index(["c0", "c1", "c2"]),
    )

    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="string"),
            ObsSchemaDescriptor(name="missing_col", dtype="string", nullable=True),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization._sanitize_obs_with_schema(adata=adata, ingest_schema=ingest_schema)

    assert list(adata.obs.columns) == ["cell_type", "missing_col"]
    assert adata.obs["missing_col"].isna().all()


def test_sanitize_obs_with_schema_preserves_column_order() -> None:
    """Verify columns are ordered according to schema."""

    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    # Create obs with columns in different order than schema
    adata.obs = pd.DataFrame(
        {"z_col": ["x", "y", "z"], "a_col": [1, 2, 3], "m_col": ["A", "B", "C"]},
        index=pd.Index(["c0", "c1", "c2"]),
    )

    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="a_col", dtype="int64"),
            ObsSchemaDescriptor(name="m_col", dtype="string"),
            ObsSchemaDescriptor(name="z_col", dtype="string"),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization._sanitize_obs_with_schema(adata=adata, ingest_schema=ingest_schema)

    # Columns should be in schema order
    assert list(adata.obs.columns) == ["a_col", "m_col", "z_col"]


# Tests for sanitize_first_adata_for_schema


@pytest.fixture()
def ingest_schema_with_full_features():
    """Create an IngestSchema with full feature set (no var columns)."""

    # Create schema var DataFrame with all features - no columns, only index
    var_df = pd.DataFrame(
        index=["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004"],
    )
    return IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="string"),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )


@pytest.fixture()
def adata_with_partial_features() -> anndata.AnnData:
    """Create an AnnData with a subset of features for schema expansion testing."""
    n_obs, n_vars = 3, 2
    X = sp.random(n_obs, n_vars, density=0.5, format="csr", dtype=np.float32)

    obs = pd.DataFrame(
        data={"cell_type": ["A", "B", "C"], "obs_id": ["obs_0", "obs_1", "obs_2"]},
        index=pd.Index(["obs_0", "obs_1", "obs_2"]),
    )
    var = pd.DataFrame(
        index=pd.Index(["ENSG0001", "ENSG0003"]),  # No columns
    )

    return anndata.AnnData(X=X, obs=obs, var=var)


def test_sanitize_first_adata_for_schema_replaces_var_from_schema(
    adata_with_partial_features: anndata.AnnData,
    ingest_schema_with_full_features,
) -> None:
    """Verify that var is replaced entirely from schema."""
    original_n_obs = adata_with_partial_features.n_obs

    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_partial_features,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Should now have all features from schema
    assert adata_with_partial_features.n_vars == 5
    assert list(adata_with_partial_features.var_names) == [
        "ENSG0000",
        "ENSG0001",
        "ENSG0002",
        "ENSG0003",
        "ENSG0004",
    ]
    # Obs count should be unchanged
    assert adata_with_partial_features.n_obs == original_n_obs
    # X shape should match
    assert adata_with_partial_features.X.shape == (original_n_obs, 5)


def test_sanitize_first_adata_for_schema_uses_schema_var_index(
    adata_with_partial_features: anndata.AnnData,
    ingest_schema_with_full_features,
) -> None:
    """Verify var index comes from schema DataFrame (no columns)."""
    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_partial_features,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Var should have no columns (only index)
    assert len(adata_with_partial_features.var.columns) == 0
    # Index should match schema
    assert list(adata_with_partial_features.var.index) == [
        "ENSG0000",
        "ENSG0001",
        "ENSG0002",
        "ENSG0003",
        "ENSG0004",
    ]


def test_sanitize_first_adata_for_schema_zero_fills_x_for_missing(
    adata_with_partial_features: anndata.AnnData,
    ingest_schema_with_full_features,
) -> None:
    """Verify X rows for missing features are zero-filled."""
    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_partial_features,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Get X as dense for easier inspection
    X_dense = (
        adata_with_partial_features.X.toarray()
        if sp.issparse(adata_with_partial_features.X)
        else adata_with_partial_features.X
    )

    # Missing features: ENSG0000, ENSG0002, ENSG0004 (indices 0, 2, 4 after reordering)
    assert np.all(X_dense[:, 0] == 0)  # ENSG0000
    assert np.all(X_dense[:, 2] == 0)  # ENSG0002
    assert np.all(X_dense[:, 4] == 0)  # ENSG0004


def test_sanitize_first_adata_for_schema_preserves_existing_data(
    adata_with_partial_features: anndata.AnnData,
    ingest_schema_with_full_features,
) -> None:
    """Verify existing X data and var values are preserved."""
    # Store original data
    original_x = adata_with_partial_features[:, ["ENSG0001", "ENSG0003"]].X.toarray().copy()

    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_partial_features,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Get X as dense
    X_dense = (
        adata_with_partial_features.X.toarray()
        if sp.issparse(adata_with_partial_features.X)
        else adata_with_partial_features.X
    )

    # Original features should be at indices 1 and 3 after reordering
    np.testing.assert_array_equal(X_dense[:, 1], original_x[:, 0])  # ENSG0001
    np.testing.assert_array_equal(X_dense[:, 3], original_x[:, 1])  # ENSG0003


def test_sanitize_first_adata_for_schema_no_op_when_all_features_present(
    ingest_schema_with_full_features,
) -> None:
    """Verify no expansion when all features are already present."""
    n_obs = 3
    X = sp.random(n_obs, 5, density=0.5, format="csr", dtype=np.float32)
    obs = pd.DataFrame(
        data={"cell_type": ["A", "B", "C"], "obs_id": ["obs_0", "obs_1", "obs_2"]},
        index=pd.Index(["obs_0", "obs_1", "obs_2"]),
    )
    var = pd.DataFrame(
        index=pd.Index(["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004"]),  # No columns
    )
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    original_shape = adata.X.shape

    sanitization.sanitize_first_adata_for_schema(
        adata=adata,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Shape should be unchanged
    assert adata.X.shape == original_shape
    assert adata.n_vars == 5
    assert len(adata.var.columns) == 0


def test_sanitize_first_adata_for_schema_calls_sanitize_for_ingest(
    adata_with_all_slots: anndata.AnnData,
) -> None:
    """Verify standard sanitization is also applied."""

    # Create schema var DataFrame matching the adata_with_all_slots fixture (no columns)
    var_df = pd.DataFrame(
        index=["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003"],
    )
    ingest_schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="string")],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_all_slots,
        ingest_schema=ingest_schema,
    )

    # Standard sanitization should have been applied
    assert len(adata_with_all_slots.obsm) == 0
    assert len(adata_with_all_slots.varm) == 0
    assert len(adata_with_all_slots.uns) == 0
    assert len(adata_with_all_slots.obsp) == 0
    assert len(adata_with_all_slots.varp) == 0
    assert len(adata_with_all_slots.layers) == 0
    assert adata_with_all_slots.raw is None
    # Var should have no columns
    assert len(adata_with_all_slots.var.columns) == 0


def test_sanitize_obs_with_schema_fills_nullable_missing_columns_numeric() -> None:
    """Verify missing nullable columns (numeric/bool) are filled with NA."""

    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    adata.obs = pd.DataFrame(
        {"cell_type": ["A", "B", "C"]},
        index=pd.Index(["c0", "c1", "c2"]),
    )

    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="string"),
            ObsSchemaDescriptor(name="int_col", dtype="int32", nullable=True),
            ObsSchemaDescriptor(name="float_col", dtype="float32", nullable=True),
            ObsSchemaDescriptor(name="bool_col", dtype="bool", nullable=True),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization._sanitize_obs_with_schema(adata=adata, ingest_schema=ingest_schema)

    # Check int column (should be Int32 or similar nullable)
    assert "int_col" in adata.obs.columns
    assert str(adata.obs["int_col"].dtype) == "Int32"
    assert adata.obs["int_col"].isna().all()

    # Check float column
    assert "float_col" in adata.obs.columns
    assert str(adata.obs["float_col"].dtype) == "float32"
    assert adata.obs["float_col"].isna().all()

    # Check bool column (should be boolean)
    assert "bool_col" in adata.obs.columns
    assert str(adata.obs["bool_col"].dtype) == "boolean"
    assert adata.obs["bool_col"].isna().all()


def test_sanitize_obs_existing_nullable_with_nans() -> None:
    """Verify existing nullable columns with NAs are safely cast to nullable types."""

    # Initial DataFrame with mixed values including NAs/None
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))
    adata.obs = pd.DataFrame(
        {
            "cell_type": ["A", "B", "C"],
            # Mix of float/nan, simulating reading from file where int column became float due to NaNs
            "total_mrna_umis": [100.0, np.nan, 200.0],
            "disease": ["d1", None, "d2"],
        },
        index=pd.Index(["c0", "c1", "c2"]),
    )

    var_df = pd.DataFrame(index=["0", "1"])
    ingest_schema = IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="string"),
            # Schema expects int32. If we cast blindly to 'int32', the NaN would fail.
            # Fix should cast this to 'Int32' instead.
            ObsSchemaDescriptor(name="total_mrna_umis", dtype="int32", nullable=True),
            ObsSchemaDescriptor(name="disease", dtype="string", nullable=True),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    sanitization.sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)

    # Verify mrna column is Int32 and preserves NA
    assert str(adata.obs["total_mrna_umis"].dtype) == "Int32"
    assert adata.obs["total_mrna_umis"].iloc[0] == 100
    assert pd.isna(adata.obs["total_mrna_umis"].iloc[1])
    assert adata.obs["total_mrna_umis"].iloc[2] == 200

    # Verify disease column is string/object and preserves NA
    # Note: 'str' dtype might be 'string' (Pandas nullable) or object depending on impl details,
    # but as long as it handles None it is fine. Our mapping uses 'string'.
    assert str(adata.obs["disease"].dtype) == "string"
    assert adata.obs["disease"].iloc[1] is pd.NA
