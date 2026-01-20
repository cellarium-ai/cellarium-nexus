"""
Sanitization functions for AnnData before SOMA ingest.

This module provides functions to remove unsupported AnnData slots
and prepare data for TileDB SOMA ingestion.
"""

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

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


def _replace_var_with_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Replace AnnData var with schema var DataFrame and reorder/expand X in-place.

    Extract feature IDs from input AnnData, validate they are in schema,
    then replace var with schema var DataFrame. Expand X to include all
    schema features (zero-fill missing), and reorder to match schema order.

    :param adata: The AnnData object to modify. Modified in-place.
    :param ingest_schema: Schema containing the full var DataFrame.
    """
    schema_var_df = ingest_schema.var_schema.to_dataframe()
    schema_features = list(schema_var_df.index)
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

    # Update adata in-place
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


def sanitize_for_ingest(*, adata: AnnData) -> None:
    """
    Sanitize AnnData for TileDB SOMA ingest in-place.

    Remove unsupported slots (obsm, varm, uns, obsp, varp, layers, raw)
    and reset index names to prepare for SOMA registration.

    :param adata: The AnnData object to sanitize.
    """
    _remove_obsm(adata=adata)
    _remove_varm(adata=adata)
    _remove_uns(adata=adata)
    _remove_obsp(adata=adata)
    _remove_varp(adata=adata)
    _remove_layers(adata=adata)
    _remove_raw(adata=adata)
    _reset_index_names(adata=adata)


def sanitize_first_adata_for_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Sanitize the first AnnData for SOMA schema creation in-place.

    Perform standard sanitization, then replace var DataFrame entirely with
    the schema var DataFrame. X matrix is expanded to include all schema
    features (zero-fill missing) and reordered to match schema feature order.

    :param adata: The AnnData object to sanitize. Modified in-place.
    :param ingest_schema: Schema containing the full var DataFrame.
    """
    sanitize_for_ingest(adata=adata)
    _replace_var_with_schema(adata=adata, ingest_schema=ingest_schema)
