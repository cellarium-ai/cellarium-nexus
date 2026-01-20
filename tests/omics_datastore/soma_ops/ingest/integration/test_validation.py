"""
Integration tests for SOMA ingest validation.

Test validate_for_ingest with real AnnData objects to verify
validation logic works end-to-end before ingestion.
"""

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import validate_for_ingest
from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaValidationError
from cellarium.nexus.shared.schemas.omics_datastore import (
    ExperimentVarSchema,
    IngestSchema,
    ObsSchemaDescriptor,
)

# Tests for validate_for_ingest


def test_validate_for_ingest_valid_anndata_passes(
    small_valid_anndata: tuple[anndata.AnnData, str],
    small_ingest_schema: IngestSchema,
) -> None:
    """Verify valid AnnData passes schema validation without error."""
    adata, _ = small_valid_anndata

    # Should not raise
    validate_for_ingest(adata=adata, schema=small_ingest_schema)


def test_validate_for_ingest_missing_required_obs_column_raises(
    small_valid_anndata: tuple[anndata.AnnData, str],
    small_ingest_schema: IngestSchema,
) -> None:
    """Verify missing required obs column raises SomaValidationError."""
    adata, _ = small_valid_anndata
    adata.obs = adata.obs.drop(columns=["cell_type"])

    with pytest.raises(SomaValidationError) as exc_info:
        validate_for_ingest(adata=adata, schema=small_ingest_schema)

    assert "cell_type" in str(exc_info.value)


def test_validate_for_ingest_unknown_var_features_raises(
    small_valid_anndata: tuple[anndata.AnnData, str],
) -> None:
    """Verify AnnData with features not in schema raises error."""
    adata, _ = small_valid_anndata

    # Create schema with different features
    var_df = pd.DataFrame(
        data={"feature_name": ["unknown_1", "unknown_2"]},
        index=["UNKNOWN_GENE_1", "UNKNOWN_GENE_2"],
    )
    schema = IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    with pytest.raises(SomaValidationError):
        validate_for_ingest(adata=adata, schema=schema)


def test_validate_for_ingest_negative_values_in_count_matrix_raises(
    tmp_path,
) -> None:
    """Verify negative values in count matrix raises SomaValidationError."""
    n_obs = 5

    # Create matrix with negative values
    X = sp.csr_matrix(np.array([[-1, 0, 1], [0, 2, 0], [0, 0, 3], [1, 0, 0], [0, 1, 0]], dtype=np.float32))

    obs = pd.DataFrame(
        data={"cell_type": ["a", "b", "c", "d", "e"]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        data={"feature_name": ["g1", "g2", "g3"]},
        index=["ENSG0001", "ENSG0002", "ENSG0003"],
    )

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    var_df = pd.DataFrame(
        data={"feature_name": ["g1", "g2", "g3"]},
        index=["ENSG0001", "ENSG0002", "ENSG0003"],
    )
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )

    with pytest.raises(SomaValidationError) as exc_info:
        validate_for_ingest(adata=adata, schema=schema)

    assert "negative" in str(exc_info.value).lower()


def test_validate_for_ingest_feature_matrix_allows_negative_values(
    tmp_path,
) -> None:
    """Verify feature_matrix validation type allows negative values."""
    n_obs = 5

    # Create matrix with negative float values
    X = sp.csr_matrix(np.array([[-1.5, 0.5, 1.2], [0, -2.3, 0], [0, 0, 3.1], [1.0, 0, 0], [0, 1.0, 0]]))

    obs = pd.DataFrame(
        data={"cell_type": ["a", "b", "c", "d", "e"]},
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        data={"feature_name": ["g1", "g2", "g3"]},
        index=["ENSG0001", "ENSG0002", "ENSG0003"],
    )

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    var_df = pd.DataFrame(
        data={"feature_name": ["g1", "g2", "g3"]},
        index=["ENSG0001", "ENSG0002", "ENSG0003"],
    )
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="feature_matrix",
    )

    # Should not raise
    validate_for_ingest(adata=adata, schema=schema)
