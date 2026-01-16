import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import sanitization


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


# Tests for sanitize_for_ingest


def test_sanitize_for_ingest_removes_all_unsupported_slots(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify all unsupported slots are removed."""
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

    sanitization.sanitize_for_ingest(adata=adata_with_all_slots)

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


def test_sanitize_for_ingest_preserves_x_obs_var(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify X, obs, and var are preserved after sanitization."""
    original_x_shape = adata_with_all_slots.X.shape
    original_obs_columns = list(adata_with_all_slots.obs.columns)
    original_var_columns = list(adata_with_all_slots.var.columns)
    original_n_obs = adata_with_all_slots.n_obs
    original_n_vars = adata_with_all_slots.n_vars

    sanitization.sanitize_for_ingest(adata=adata_with_all_slots)

    assert adata_with_all_slots.X.shape == original_x_shape
    assert list(adata_with_all_slots.obs.columns) == original_obs_columns
    assert list(adata_with_all_slots.var.columns) == original_var_columns
    assert adata_with_all_slots.n_obs == original_n_obs
    assert adata_with_all_slots.n_vars == original_n_vars


def test_sanitize_for_ingest_operates_in_place(adata_with_all_slots: anndata.AnnData) -> None:
    """Verify sanitization operates in-place on the same object."""
    original_id = id(adata_with_all_slots)

    sanitization.sanitize_for_ingest(adata=adata_with_all_slots)

    assert id(adata_with_all_slots) == original_id


def test_sanitize_for_ingest_handles_minimal_adata() -> None:
    """Verify sanitization handles minimal AnnData with only X."""
    adata = anndata.AnnData(X=sp.csr_matrix((3, 2)))

    sanitization.sanitize_for_ingest(adata=adata)

    assert len(adata.obsm) == 0
    assert len(adata.varm) == 0
    assert len(adata.uns) == 0
    assert len(adata.obsp) == 0
    assert len(adata.varp) == 0
    assert len(adata.layers) == 0
    assert adata.raw is None


# Tests for sanitize_first_adata_for_schema


@pytest.fixture()
def ingest_schema_with_full_features():
    """Create an IngestSchema with full feature set and var columns."""
    from cellarium.nexus.shared.schemas.omics_datastore import (
        ExperimentVarFeatures,
        IngestSchema,
        ObsSchemaDescriptor,
        VarSchemaDescriptor,
    )

    return IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="str"),
        ],
        var_columns=[
            VarSchemaDescriptor(name="gene_name", dtype="str"),
            VarSchemaDescriptor(name="feature_length", dtype="int"),
        ],
        var_features=ExperimentVarFeatures(
            features=["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004"],
            is_subset=True,
        ),
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
        data={"gene_name": ["gene_0", "gene_1"], "var_id": ["ENSG0001", "ENSG0003"]},
        index=pd.Index(["ENSG0001", "ENSG0003"]),
    )

    return anndata.AnnData(X=X, obs=obs, var=var)


def test_sanitize_first_adata_for_schema_expands_features(
    adata_with_partial_features: anndata.AnnData,
    ingest_schema_with_full_features,
) -> None:
    """Verify that missing features are added to var and X."""
    original_n_obs = adata_with_partial_features.n_obs

    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_partial_features,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Should now have all features
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


def test_sanitize_first_adata_for_schema_adds_var_columns(
    adata_with_partial_features: anndata.AnnData,
    ingest_schema_with_full_features,
) -> None:
    """Verify var columns are added with default value 0."""
    sanitization.sanitize_first_adata_for_schema(
        adata=adata_with_partial_features,
        ingest_schema=ingest_schema_with_full_features,
    )

    # Both var columns should be present
    assert "gene_name" in adata_with_partial_features.var.columns
    assert "feature_length" in adata_with_partial_features.var.columns

    # Missing features should have 0 for the new var columns
    missing_feature_row = adata_with_partial_features.var.loc["ENSG0000"]
    assert missing_feature_row["gene_name"] == 0
    assert missing_feature_row["feature_length"] == 0


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
        data={
            "gene_name": ["g0", "g1", "g2", "g3", "g4"],
            "var_id": ["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004"],
        },
        index=pd.Index(["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004"]),
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


def test_sanitize_first_adata_for_schema_calls_sanitize_for_ingest(
    adata_with_all_slots: anndata.AnnData,
) -> None:
    """Verify standard sanitization is also applied."""
    from cellarium.nexus.shared.schemas.omics_datastore import (
        ExperimentVarFeatures,
        IngestSchema,
        ObsSchemaDescriptor,
        VarSchemaDescriptor,
    )

    # Use features that match the adata_with_all_slots fixture
    ingest_schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str")],
        var_columns=[VarSchemaDescriptor(name="gene_name", dtype="str")],
        var_features=ExperimentVarFeatures(
            features=["ENSG0000", "ENSG0001", "ENSG0002", "ENSG0003"],
            is_subset=True,
        ),
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
