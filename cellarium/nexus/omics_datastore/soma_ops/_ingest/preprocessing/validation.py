"""
Validation functions for AnnData before SOMA ingest.

This module provides validation for obs schema (dtype castability, nullability),
var schema (feature membership), and X matrix (count vs feature matrix).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing.constants import PANDAS_NULLABLE_DTYPES
from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaValidationError
from cellarium.nexus.shared.schemas.omics_datastore import ExperimentVarSchema, IngestSchema, ObsSchemaDescriptor

# Type alias for X matrix
XMatrix = sp.spmatrix | np.ndarray[Any, Any]
XMatrixOptional = XMatrix | None


def _try_cast_column(*, series: pd.Series[Any], target_dtype: Any) -> tuple[bool, str | None]:
    """
    Attempt to cast a pandas Series to the target dtype.

    :param series: The pandas Series to cast.
    :param target_dtype: The target dtype as a string.

    :returns: Tuple of (success, error_message). If success is True, error_message is None.
    """
    try:
        # Use nullable dtype for compatibility check if available (e.g. "int32" -> "Int32")
        # so that NaNs in nullable columns don't cause cast failure.
        cast_dtype = PANDAS_NULLABLE_DTYPES.get(target_dtype, target_dtype)
        series.astype(cast_dtype)  # type: ignore[arg-type]
        return True, None
    except (ValueError, TypeError) as e:
        return False, str(e)


def _validate_obs_column_presence(
    *,
    obs_df: pd.DataFrame,
    obs_columns: list[ObsSchemaDescriptor],
) -> list[str]:
    """
    Validate that required (non-nullable) columns exist in obs.

    :param obs_df: The obs DataFrame from AnnData.
    :param obs_columns: List of obs descriptors to validate against.

    :returns: List of error messages for missing non-nullable columns.
    """
    errors: list[str] = []
    for col_schema in obs_columns:
        if col_schema.name not in obs_df.columns and not col_schema.nullable:
            errors.append(f"obs: Required column '{col_schema.name}' is missing")
    return errors


def _validate_obs_column_dtypes(
    *,
    obs_df: pd.DataFrame,
    obs_columns: list[ObsSchemaDescriptor],
) -> list[str]:
    """
    Validate that existing obs columns can be cast to their target dtypes.

    :param obs_df: The obs DataFrame from AnnData.
    :param obs_columns: List of obs descriptors to validate against.

    :returns: List of error messages for columns that cannot be cast.
    """
    errors: list[str] = []
    for col_schema in obs_columns:
        # Only validate dtype for columns that are present
        # Missing nullable columns are allowed (handled in presence validation)
        if col_schema.name in obs_df.columns:
            success, error_msg = _try_cast_column(series=obs_df[col_schema.name], target_dtype=col_schema.dtype)
            if not success:
                errors.append(f"obs: Column '{col_schema.name}' cannot be cast to '{col_schema.dtype}': {error_msg}")
    return errors


def validate_obs(*, obs_df: pd.DataFrame, obs_columns: list[ObsSchemaDescriptor]) -> list[str]:
    """
    Validate obs DataFrame against schema.

    Check column presence (non-nullable columns must exist) and dtype castability
    (existing columns must be castable to target dtype).

    :param obs_df: The obs DataFrame from AnnData.
    :param obs_columns: List of obs descriptors to validate against.

    :returns: List of all validation error messages.
    """
    errors: list[str] = []
    errors.extend(_validate_obs_column_presence(obs_df=obs_df, obs_columns=obs_columns))
    errors.extend(_validate_obs_column_dtypes(obs_df=obs_df, obs_columns=obs_columns))
    return errors


def validate_var(*, var_df: pd.DataFrame, var_schema: ExperimentVarSchema) -> list[str]:
    """
    Validate var DataFrame against schema.

    Check that var index is string type and all features in var index are
    present in the schema. Extra features not in schema are always rejected.
    If is_subset=False, also verify that all schema features are present.

    :param var_df: The var DataFrame from AnnData.
    :param var_schema: The experiment var schema to validate against.

    :returns: List of all validation error messages.
    """
    errors: list[str] = []

    # Check var index is string type
    if not pd.api.types.is_string_dtype(var_df.index) and not pd.api.types.is_object_dtype(var_df.index):
        errors.append(f"var: Index must be string type, got '{var_df.index.dtype}'")

    # Get input features from var index
    input_features = set(var_df.index)
    schema_features = set(var_schema.get_feature_ids())

    # Check all input features are in schema (reject extra features always)
    unknown_features = input_features - schema_features
    if unknown_features:
        sample = list(unknown_features)[:5]
        errors.append(f"var: Found {len(unknown_features)} feature(s) not in schema. " f"Examples: {sample}")

    # If is_subset=False, all schema features must be present
    if not var_schema.is_subset:
        missing_features = schema_features - input_features
        if missing_features:
            sample = list(missing_features)[:5]
            errors.append(
                f"var: Missing {len(missing_features)} required feature(s) from schema. " f"Examples: {sample}"
            )

    return errors


def _get_x_values(*, X: XMatrix) -> np.ndarray[Any, Any]:
    """
    Extract values array from X matrix.

    :param X: The X matrix (sparse or dense).

    :returns: 1D numpy array of values.
    """
    if sp.issparse(X):
        return np.asarray(X.data)  # type: ignore[union-attr]
    return np.asarray(X).ravel()


def _validate_x_finiteness(*, X: XMatrix) -> list[str]:
    """
    Validate that X contains only finite values (no NaN, inf, -inf).

    :param X: The X matrix (sparse or dense).

    :returns: List of error messages.
    """
    errors: list[str] = []
    values = _get_x_values(X=X)

    if not np.all(np.isfinite(values)):
        nan_count = np.sum(np.isnan(values))
        inf_count = np.sum(np.isinf(values))
        parts = []
        if nan_count > 0:
            parts.append(f"{nan_count} NaN")
        if inf_count > 0:
            parts.append(f"{inf_count} infinite")
        errors.append(f"X: Contains non-finite values: {', '.join(parts)}")

    return errors


def _validate_x_shape(*, X: XMatrix, n_obs: int, n_vars: int) -> list[str]:
    """
    Validate that X shape matches expected dimensions.

    :param X: The X matrix (sparse or dense).
    :param n_obs: Expected number of observations (rows).
    :param n_vars: Expected number of variables (columns).

    :returns: List of error messages.
    """
    errors: list[str] = []
    if X.shape != (n_obs, n_vars):
        errors.append(f"X: Shape {X.shape} does not match expected ({n_obs}, {n_vars})")
    return errors


def _validate_x_non_negativity(*, X: XMatrix) -> list[str]:
    """
    Validate that X contains only non-negative values.

    :param X: The X matrix (sparse or dense).

    :returns: List of error messages.
    """
    errors: list[str] = []
    values = _get_x_values(X=X)

    if np.any(values < 0):
        neg_count = np.sum(values < 0)
        errors.append(f"X: Contains {neg_count} negative value(s)")

    return errors


def _validate_x_integer_like(*, X: XMatrix) -> list[str]:
    """
    Validate that X contains integer-like values.

    Allow floating-point representation of integers (e.g., 3.0 is valid).

    :param X: The X matrix (sparse or dense).

    :returns: List of error messages.
    """
    errors: list[str] = []
    values = _get_x_values(X=X)

    # Skip check for integer dtypes
    if np.issubdtype(values.dtype, np.integer):
        return errors

    # Check if float values are integer-like
    if not np.allclose(values, np.round(values)):
        non_int_mask = ~np.isclose(values, np.round(values))
        non_int_count = np.sum(non_int_mask)
        sample_values = values[non_int_mask][:5]
        errors.append(f"X: Contains {non_int_count} non-integer value(s). " f"Examples: {sample_values.tolist()}")

    return errors


def _validate_x_int_range(*, X: XMatrix) -> list[str]:
    """
    Validate that X values are within int32 range.

    :param X: The X matrix (sparse or dense).

    :returns: List of error messages.
    """
    errors: list[str] = []
    values = _get_x_values(X=X)

    int32_min = np.iinfo(np.int32).min
    int32_max = np.iinfo(np.int32).max

    if np.any(values < int32_min) or np.any(values > int32_max):
        out_of_range = np.sum((values < int32_min) | (values > int32_max))
        errors.append(f"X: Contains {out_of_range} value(s) outside int32 range " f"[{int32_min}, {int32_max}]")

    return errors


def _validate_count_matrix(*, X: XMatrix) -> list[str]:
    """
    Validate X as a count matrix (raw counts).

    Check finiteness, non-negativity, integer-like values, and int32 range.

    :param X: The X matrix (sparse or dense).

    :returns: List of error messages.
    """
    errors: list[str] = []
    errors.extend(_validate_x_finiteness(X=X))
    errors.extend(_validate_x_non_negativity(X=X))
    errors.extend(_validate_x_integer_like(X=X))
    errors.extend(_validate_x_int_range(X=X))
    return errors


def _validate_feature_matrix(*, X: XMatrix) -> list[str]:
    """
    Validate X as a feature matrix (allows any numeric values).

    Check finiteness only.

    :param X: The X matrix (sparse or dense).

    :returns: List of error messages.
    """
    errors: list[str] = []
    errors.extend(_validate_x_finiteness(X=X))
    return errors


def validate_x(
    *,
    X: XMatrixOptional,
    validation_type: str,
    n_obs: int,
    n_vars: int,
) -> list[str]:
    """
    Validate X matrix against schema.

    Always validate shape, then dispatch to count_matrix or feature_matrix validation.

    :param X: The X matrix (sparse or dense).
    :param validation_type: Either "count_matrix" or "feature_matrix".
    :param n_obs: Expected number of observations (rows).
    :param n_vars: Expected number of variables (columns).

    :returns: List of all validation error messages.
    """
    errors: list[str] = []

    if X is None:
        errors.append("X: Matrix is None")
        return errors

    errors.extend(_validate_x_shape(X=X, n_obs=n_obs, n_vars=n_vars))

    if validation_type == "count_matrix":
        errors.extend(_validate_count_matrix(X=X))
    elif validation_type == "feature_matrix":
        errors.extend(_validate_feature_matrix(X=X))
    else:
        errors.append(f"X: Unknown validation type '{validation_type}'")

    return errors


def validate_for_ingest(*, adata: AnnData, schema: IngestSchema) -> None:
    """
    Validate AnnData for TileDB SOMA ingest.

    Run all validations (obs, var, X) and collect all errors. Raise SomaValidationError
    with full error list if any validation fails.

    :param adata: The AnnData object to validate.
    :param schema: The ingest schema to validate against.

    :raises SomaValidationError: If any validation check fails.
    """
    errors = []

    errors.extend(validate_obs(obs_df=adata.obs, obs_columns=schema.obs_columns))
    errors.extend(validate_var(var_df=adata.var, var_schema=schema.var_schema))
    errors.extend(
        validate_x(X=adata.X, validation_type=schema.x_validation_type, n_obs=adata.n_obs, n_vars=adata.n_vars)
    )

    if errors:
        raise SomaValidationError(errors=errors)
