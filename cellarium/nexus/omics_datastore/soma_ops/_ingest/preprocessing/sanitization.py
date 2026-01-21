"""
Sanitization functions for AnnData before SOMA ingest.

This module provides functions to remove unsupported AnnData slots
and prepare data for TileDB SOMA ingestion.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing.constants import PANDAS_NULLABLE_DTYPES
from cellarium.nexus.shared.schemas.omics_datastore import IngestSchema


def _remove_obsm(*, adata: AnnData) -> None:
    """
    Remove all obsm (observation embeddings) from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    for key in list(adata.obsm.keys()):
        del adata.obsm[key]


def _remove_varm(*, adata: AnnData) -> None:
    """
    Remove all varm (variable embeddings) from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    for key in list(adata.varm.keys()):
        del adata.varm[key]


def _remove_uns(*, adata: AnnData) -> None:
    """
    Remove all uns (unstructured annotations) from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    adata.uns.clear()


def _remove_obsp(*, adata: AnnData) -> None:
    """
    Remove all obsp (observation pairwise data) from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    for key in list(adata.obsp.keys()):
        del adata.obsp[key]


def _remove_varp(*, adata: AnnData) -> None:
    """
    Remove all varp (variable pairwise data) from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    for key in list(adata.varp.keys()):
        del adata.varp[key]


def _remove_layers(*, adata: AnnData) -> None:
    """
    Remove all layers from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    for key in list(adata.layers.keys()):
        del adata.layers[key]


def _remove_raw(*, adata: AnnData) -> None:
    """
    Remove raw data from AnnData in-place.

    :param adata: The AnnData object to sanitize.
    """
    if adata.raw is not None:
        del adata.raw


def _reset_index_names(*, adata: AnnData) -> None:
    """
    Reset obs and var index names to None.

    SOMA uses its own indexing scheme (soma_joinid), so original
    index names are discarded.

    :param adata: The AnnData object to sanitize.
    """
    adata.obs.index.name = None
    adata.var.index.name = None


def _sanitize_var_metadata_with_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Sanitize var DataFrame metadata according to schema in-place.

    Replace var DataFrame with the subset of schema features that match
    the AnnData's features. This ensures var columns and their types
    match the SOMA experiment schema for the existing features.

    :param adata: The AnnData object to sanitize.
    :param ingest_schema: Schema containing the full var DataFrame.
    """
    schema_var_df = ingest_schema.var_schema.to_dataframe()
    # Align schema var to adata features (preserving adata order)
    # Validation ensures all adata features exist in schema
    adata.var = schema_var_df.reindex(adata.var.index)


def _sanitize_obs_with_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Sanitize obs DataFrame according to schema in-place.

    Filter obs columns to only those defined in schema and cast each column
    to its target dtype. Columns not in schema are dropped. Missing nullable
    columns are filled with NA values.

    :param adata: The AnnData object to sanitize.
    :param ingest_schema: Schema containing obs column descriptors.
    """
    obs_columns = ingest_schema.obs_columns
    schema_column_names = [col.name for col in obs_columns]

    # Build new obs DataFrame with only schema columns
    new_obs_data: dict[str, pd.Series] = {}
    for col_schema in obs_columns:
        col_name = col_schema.name
        if col_name in adata.obs.columns:
            # Cast to target dtype
            if col_schema.nullable:
                dtype = PANDAS_NULLABLE_DTYPES.get(col_schema.dtype, col_schema.dtype)
                new_obs_data[col_name] = adata.obs[col_name].astype(dtype)
            else:
                new_obs_data[col_name] = adata.obs[col_name].astype(col_schema.dtype)
        elif col_schema.nullable:
            # Create column with NA values for nullable missing columns
            dtype = PANDAS_NULLABLE_DTYPES.get(col_schema.dtype, col_schema.dtype)
            new_obs_data[col_name] = pd.Series(
                [pd.NA] * adata.n_obs,
                index=adata.obs.index,
                dtype=dtype,
            )
        # Non-nullable missing columns are skipped - validation should catch this

    # Replace obs with sanitized DataFrame, preserving index
    new_obs = pd.DataFrame(new_obs_data, index=adata.obs.index)
    # Ensure column order matches schema order
    new_obs = new_obs[[col for col in schema_column_names if col in new_obs.columns]]
    adata.obs = new_obs


def _replace_var_with_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Replace AnnData var with schema var DataFrame and reorder/expand X in-place.

    Extract feature IDs from input AnnData, then replace var with the full
    schema var DataFrame (including all columns). Expand X to include all
    schema features (zero-fill missing), and reorder to match schema order.

    :param adata: The AnnData object to modify. Modified in-place.
    :param ingest_schema: Schema containing the full var DataFrame.
    """
    schema_features = ingest_schema.var_schema.get_feature_ids()
    input_features = list(adata.var_names)

    # Build mapping from schema feature to column index
    schema_feature_to_idx = {f: i for i, f in enumerate(schema_features)}

    n_cells = adata.n_obs
    n_schema_features = len(schema_features)

    # Determine X dtype
    x_dtype = adata.X.dtype if adata.X is not None else np.float32  # type: ignore[union-attr]

    # Create new X matrix with schema dimensions
    if sp.issparse(adata.X):
        # Build COO data for new sparse matrix
        old_X = adata.X.tocoo()  # type: ignore[union-attr]

        # Map old column indices to new column indices
        old_to_new_col = {}
        for old_idx, feature in enumerate(input_features):
            if feature in schema_feature_to_idx:
                old_to_new_col[old_idx] = schema_feature_to_idx[feature]

        # Filter and remap
        new_rows = []
        new_cols = []
        new_data = []
        for i, (row, col, val) in enumerate(zip(old_X.row, old_X.col, old_X.data)):
            if col in old_to_new_col:
                new_rows.append(row)
                new_cols.append(old_to_new_col[col])
                new_data.append(val)

        new_X = sp.csr_matrix(
            (new_data, (new_rows, new_cols)),
            shape=(n_cells, n_schema_features),
            dtype=x_dtype,
        )
    else:
        # Dense matrix
        new_X = np.zeros((n_cells, n_schema_features), dtype=x_dtype)
        for old_idx, feature in enumerate(input_features):
            if feature in schema_feature_to_idx:
                new_idx = schema_feature_to_idx[feature]
                new_X[:, new_idx] = adata.X[:, old_idx]  # type: ignore[index]

    # Update adata in-place with var from schema (includes all columns)
    schema_var_df = ingest_schema.var_schema.to_dataframe()
    adata._init_as_actual(
        X=new_X,
        obs=adata.obs,
        var=schema_var_df,
        uns={},
        obsm={},
        varm={},
        obsp={},
        varp={},
        layers={},
        raw=None,
    )


def sanitize_for_ingest(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Sanitize AnnData for TileDB SOMA ingest in-place.

    Remove unsupported slots (obsm, varm, uns, obsp, varp, layers, raw),
    sanitize obs/var columns according to schema, and reset index names
    to prepare for SOMA ingestion.

    :param adata: The AnnData object to sanitize.
    :param ingest_schema: Schema containing obs column descriptors for filtering and casting.
    """
    _remove_obsm(adata=adata)
    _remove_varm(adata=adata)
    _remove_uns(adata=adata)
    _remove_obsp(adata=adata)
    _remove_varp(adata=adata)
    _remove_layers(adata=adata)
    _remove_raw(adata=adata)
    _sanitize_var_metadata_with_schema(adata=adata, ingest_schema=ingest_schema)
    _sanitize_obs_with_schema(adata=adata, ingest_schema=ingest_schema)
    _reset_index_names(adata=adata)


def sanitize_first_adata_for_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Sanitize the first AnnData for SOMA schema creation in-place.

    Perform standard sanitization (including obs sanitization according to schema),
    then replace var DataFrame entirely with the schema var DataFrame. X matrix is
    expanded to include all schema features (zero-fill missing) and reordered to
    match schema feature order.

    :param adata: The AnnData object to sanitize. Modified in-place.
    :param ingest_schema: Schema containing the full var DataFrame and obs column descriptors.
    """
    sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)
    _replace_var_with_schema(adata=adata, ingest_schema=ingest_schema)
