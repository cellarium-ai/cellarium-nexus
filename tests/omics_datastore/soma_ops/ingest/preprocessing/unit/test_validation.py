import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import validation
from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaValidationError
from cellarium.nexus.shared.schemas.omics_datastore import ExperimentVarFeatures, IngestSchema, ObsSchemaDescriptor

# Fixtures


@pytest.fixture()
def obs_columns() -> list[ObsSchemaDescriptor]:
    """
    Create a list of obs column descriptors for testing.

    :return: List of ObsDescriptor objects.
    """
    return [
        ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False),
        ObsSchemaDescriptor(name="sample_id", dtype="str", nullable=False),
        ObsSchemaDescriptor(name="age", dtype="int64", nullable=True),
        ObsSchemaDescriptor(name="score", dtype="float64", nullable=True),
    ]


@pytest.fixture()
def var_schema() -> ExperimentVarFeatures:
    """
    Create a var schema for testing.

    :return: ExperimentVarSchema object.
    """
    return ExperimentVarFeatures(
        features=["ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004", "ENSG0005"],
        is_subset=True,
    )


@pytest.fixture()
def valid_obs_df() -> pd.DataFrame:
    """
    Create a valid obs DataFrame for testing.

    :return: Valid obs DataFrame.
    """
    return pd.DataFrame(
        data={
            "cell_type": ["A", "B", "C"],
            "sample_id": ["S1", "S2", "S3"],
            "age": [25, 30, 35],
            "score": [0.5, 0.8, 0.9],
        },
        index=["cell_0", "cell_1", "cell_2"],
    )


@pytest.fixture()
def valid_var_df() -> pd.DataFrame:
    """
    Create a valid var DataFrame for testing.

    :return: Valid var DataFrame.
    """
    return pd.DataFrame(
        data={"gene_name": ["gene_1", "gene_2", "gene_3"]},
        index=["ENSG0001", "ENSG0002", "ENSG0003"],
    )


@pytest.fixture()
def valid_x_count_matrix() -> sp.csr_matrix:
    """
    Create a valid count matrix for testing.

    :return: Valid CSR count matrix.
    """
    data = np.array([1, 2, 3, 4, 5], dtype=np.int32)
    rows = np.array([0, 0, 1, 2, 2])
    cols = np.array([0, 1, 2, 0, 2])
    return sp.csr_matrix((data, (rows, cols)), shape=(3, 3))


# Tests for validate_obs


def test_validate_obs_valid_dataframe(valid_obs_df: pd.DataFrame, obs_columns: list[ObsSchemaDescriptor]) -> None:
    """Verify valid obs DataFrame passes validation."""
    errors = validation.validate_obs(obs_df=valid_obs_df, obs_columns=obs_columns)

    assert errors == []


def test_validate_obs_missing_required_column(
    valid_obs_df: pd.DataFrame, obs_columns: list[ObsSchemaDescriptor]
) -> None:
    """Verify missing required column is detected."""
    obs_df = valid_obs_df.drop(columns=["cell_type"])

    errors = validation.validate_obs(obs_df=obs_df, obs_columns=obs_columns)

    assert len(errors) == 1
    assert "cell_type" in errors[0]
    assert "missing" in errors[0].lower()


def test_validate_obs_missing_nullable_column(
    valid_obs_df: pd.DataFrame, obs_columns: list[ObsSchemaDescriptor]
) -> None:
    """Verify missing nullable column is allowed."""
    obs_df = valid_obs_df.drop(columns=["age"])

    errors = validation.validate_obs(obs_df=obs_df, obs_columns=obs_columns)

    assert errors == []


def test_validate_obs_invalid_dtype_cast(valid_obs_df: pd.DataFrame, obs_columns: list[ObsSchemaDescriptor]) -> None:
    """Verify column that cannot be cast to target dtype is detected."""
    obs_df = valid_obs_df.copy()
    obs_df["age"] = ["young", "middle", "old"]  # String values cannot cast to int64

    errors = validation.validate_obs(obs_df=obs_df, obs_columns=obs_columns)

    assert len(errors) == 1
    assert "age" in errors[0]
    assert "cast" in errors[0].lower()


def test_validate_obs_multiple_errors(obs_columns: list[ObsSchemaDescriptor]) -> None:
    """Verify multiple errors are collected."""
    obs_df = pd.DataFrame(
        data={
            "sample_id": ["S1", "S2"],
            "age": ["young", "old"],  # Cannot cast to int64
        },
        index=["cell_0", "cell_1"],
    )

    errors = validation.validate_obs(obs_df=obs_df, obs_columns=obs_columns)

    assert len(errors) == 2  # Missing cell_type + invalid age dtype


# Tests for validate_var


def test_validate_var_valid_dataframe(valid_var_df: pd.DataFrame, var_schema: ExperimentVarFeatures) -> None:
    """Verify valid var DataFrame passes validation."""
    errors = validation.validate_var(var_df=valid_var_df, var_schema=var_schema)

    assert errors == []


def test_validate_var_non_string_index(var_schema: ExperimentVarFeatures) -> None:
    """Verify non-string var index is detected."""
    var_df = pd.DataFrame(
        data={"gene_name": ["gene_1", "gene_2"]},
        index=[1, 2],  # Integer index
    )

    errors = validation.validate_var(var_df=var_df, var_schema=var_schema)

    assert len(errors) >= 1
    assert any("string" in e.lower() for e in errors)


def test_validate_var_unknown_features(var_schema: ExperimentVarFeatures) -> None:
    """Verify features not in schema are detected."""
    var_df = pd.DataFrame(
        data={"gene_name": ["gene_1", "gene_2"]},
        index=["ENSG9999", "ENSG8888"],  # Not in schema
    )

    errors = validation.validate_var(var_df=var_df, var_schema=var_schema)

    assert len(errors) == 1
    assert "not in schema" in errors[0].lower()


def test_validate_var_subset_allowed(var_schema: ExperimentVarFeatures) -> None:
    """Verify subset of features is allowed when is_subset=True."""
    var_df = pd.DataFrame(
        data={"gene_name": ["gene_1"]},
        index=["ENSG0001"],  # Only one of five features
    )

    errors = validation.validate_var(var_df=var_df, var_schema=var_schema)

    assert errors == []


def test_validate_var_exact_match_required() -> None:
    """Verify exact match is required when is_subset=False."""
    var_schema = ExperimentVarFeatures(
        features=["ENSG0001", "ENSG0002", "ENSG0003"],
        is_subset=False,
    )
    var_df = pd.DataFrame(
        data={"gene_name": ["gene_1"]},
        index=["ENSG0001"],  # Missing ENSG0002, ENSG0003
    )

    errors = validation.validate_var(var_df=var_df, var_schema=var_schema)

    assert len(errors) == 1
    assert "missing" in errors[0].lower()


# Tests for validate_x


def test_validate_x_valid_count_matrix(valid_x_count_matrix: sp.csr_matrix) -> None:
    """Verify valid count matrix passes validation."""
    errors = validation.validate_x(
        X=valid_x_count_matrix,
        validation_type="count_matrix",
        n_obs=3,
        n_vars=3,
    )

    assert errors == []


def test_validate_x_none_matrix() -> None:
    """Verify None matrix is detected."""
    errors = validation.validate_x(
        X=None,
        validation_type="count_matrix",
        n_obs=3,
        n_vars=3,
    )

    assert len(errors) == 1
    assert "None" in errors[0]


def test_validate_x_wrong_shape(valid_x_count_matrix: sp.csr_matrix) -> None:
    """Verify wrong shape is detected."""
    errors = validation.validate_x(
        X=valid_x_count_matrix,
        validation_type="count_matrix",
        n_obs=5,  # Expected 5 but matrix is 3x3
        n_vars=5,
    )

    assert len(errors) >= 1
    assert any("shape" in e.lower() for e in errors)


def test_validate_x_nan_values() -> None:
    """Verify NaN values are detected."""
    X = sp.csr_matrix(np.array([[1.0, np.nan], [3.0, 4.0]]))

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert any("NaN" in e for e in errors)


def test_validate_x_infinite_values() -> None:
    """Verify infinite values are detected."""
    X = sp.csr_matrix(np.array([[1.0, np.inf], [3.0, 4.0]]))

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert any("infinite" in e for e in errors)


def test_validate_x_negative_values_count_matrix() -> None:
    """Verify negative values are detected for count matrix."""
    X = sp.csr_matrix(np.array([[1, -2], [3, 4]]))

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert any("negative" in e.lower() for e in errors)


def test_validate_x_non_integer_values_count_matrix() -> None:
    """Verify non-integer values are detected for count matrix."""
    X = sp.csr_matrix(np.array([[1.5, 2.0], [3.0, 4.0]]))

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert any("non-integer" in e.lower() for e in errors)


def test_validate_x_integer_like_floats_allowed() -> None:
    """Verify integer-like float values are allowed for count matrix."""
    X = sp.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]))

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert not any("non-integer" in e.lower() for e in errors)


def test_validate_x_out_of_int32_range() -> None:
    """Verify values outside int32 range are detected."""
    large_value = np.iinfo(np.int32).max + 1
    X = sp.csr_matrix(np.array([[large_value, 1], [2, 3]], dtype=np.int64))

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert any("int32" in e.lower() for e in errors)


def test_validate_x_feature_matrix_allows_negative() -> None:
    """Verify feature matrix allows negative values."""
    X = sp.csr_matrix(np.array([[-1.5, 2.0], [3.0, -4.0]]))

    errors = validation.validate_x(
        X=X,
        validation_type="feature_matrix",
        n_obs=2,
        n_vars=2,
    )

    # Should only check finiteness for feature matrix
    assert not any("negative" in e.lower() for e in errors)


def test_validate_x_feature_matrix_allows_floats() -> None:
    """Verify feature matrix allows non-integer float values."""
    X = sp.csr_matrix(np.array([[1.5, 2.7], [3.14, 4.99]]))

    errors = validation.validate_x(
        X=X,
        validation_type="feature_matrix",
        n_obs=2,
        n_vars=2,
    )

    assert errors == []


def test_validate_x_unknown_validation_type() -> None:
    """Verify unknown validation type is detected."""
    X = sp.csr_matrix(np.array([[1, 2], [3, 4]]))

    errors = validation.validate_x(
        X=X,
        validation_type="unknown_type",
        n_obs=2,
        n_vars=2,
    )

    assert any("unknown" in e.lower() for e in errors)


def test_validate_x_dense_matrix() -> None:
    """Verify dense matrix is handled correctly."""
    X = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int32)

    errors = validation.validate_x(
        X=X,
        validation_type="count_matrix",
        n_obs=2,
        n_vars=3,
    )

    assert errors == []


# Tests for validate_for_ingest


def test_validate_for_ingest_valid_adata(
    valid_obs_df: pd.DataFrame,
    valid_var_df: pd.DataFrame,
    valid_x_count_matrix: sp.csr_matrix,
    obs_columns: list[ObsSchemaDescriptor],
    var_schema: ExperimentVarFeatures,
) -> None:
    """Verify valid AnnData passes validation without raising."""
    adata = anndata.AnnData(X=valid_x_count_matrix, obs=valid_obs_df, var=valid_var_df)
    schema = IngestSchema(
        obs_columns=obs_columns,
        var_columns=[],
        var_features=var_schema,
        x_validation_type="count_matrix",
    )

    # Should not raise
    validation.validate_for_ingest(adata=adata, schema=schema)


def test_validate_for_ingest_raises_on_errors(
    valid_var_df: pd.DataFrame,
    valid_x_count_matrix: sp.csr_matrix,
    obs_columns: list[ObsSchemaDescriptor],
    var_schema: ExperimentVarFeatures,
) -> None:
    """Verify SomaValidationError is raised with all errors."""
    # Create invalid obs (missing required column)
    obs_df = pd.DataFrame(
        data={"sample_id": ["S1", "S2", "S3"]},
        index=["cell_0", "cell_1", "cell_2"],
    )

    adata = anndata.AnnData(X=valid_x_count_matrix, obs=obs_df, var=valid_var_df)
    schema = IngestSchema(
        obs_columns=obs_columns,
        var_columns=[],
        var_features=var_schema,
        x_validation_type="count_matrix",
    )

    with pytest.raises(SomaValidationError) as exc_info:
        validation.validate_for_ingest(adata=adata, schema=schema)

    assert len(exc_info.value.errors) >= 1
    assert "cell_type" in str(exc_info.value)


def test_validate_for_ingest_collects_all_errors(
    obs_columns: list[ObsSchemaDescriptor],
) -> None:
    """Verify all errors from obs, var, and X are collected."""
    # Create AnnData with multiple validation issues
    obs_df = pd.DataFrame(
        data={"sample_id": ["S1", "S2"]},  # Missing cell_type
        index=["cell_0", "cell_1"],
    )
    var_df = pd.DataFrame(
        data={"gene_name": ["gene_1", "gene_2"]},
        index=["UNKNOWN1", "UNKNOWN2"],  # Not in schema
    )
    X = sp.csr_matrix(np.array([[-1, np.nan], [3, 4]]))  # Negative + NaN

    adata = anndata.AnnData(X=X, obs=obs_df, var=var_df)
    schema = IngestSchema(
        obs_columns=obs_columns,
        var_columns=[],
        var_features=ExperimentVarFeatures(features=["ENSG0001", "ENSG0002"], is_subset=True),
        x_validation_type="count_matrix",
    )

    with pytest.raises(SomaValidationError) as exc_info:
        validation.validate_for_ingest(adata=adata, schema=schema)

    # Should have errors from obs, var, and X
    errors = exc_info.value.errors
    assert any("cell_type" in e for e in errors)  # obs error
    assert any("not in schema" in e.lower() for e in errors)  # var error
    assert any("NaN" in e or "negative" in e.lower() for e in errors)  # X errors
