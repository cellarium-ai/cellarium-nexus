"""
Sanitization functions for AnnData before SOMA ingest.

This module provides functions to remove unsupported AnnData slots
and prepare data for TileDB SOMA ingestion.
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData, concat

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


def _expand_to_full_feature_schema(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Expand AnnData var to include all features from the ingest schema in-place.

    Add missing features with zero-filled X rows and var columns defaulted to 0.
    Reorder features to match the schema feature order.

    :param adata: The AnnData object to expand. Modified in-place.
    :param ingest_schema: Schema containing the full feature set and var columns.
    """
    all_features = ingest_schema.var_features.features
    current_features = set(adata.var_names)
    missing_features = [f for f in all_features if f not in current_features]

    # Ensure existing var has all required columns (default to 0 if missing)
    for var_col in ingest_schema.var_columns:
        if var_col.name not in adata.var.columns:
            adata.var[var_col.name] = 0

    if missing_features:
        n_cells = adata.n_obs
        n_missing = len(missing_features)

        # Create var DataFrame for missing features with columns defaulted to 0
        missing_var_data: dict[str, list] = {"var_id": missing_features}
        for var_col in ingest_schema.var_columns:
            missing_var_data[var_col.name] = [0] * n_missing

        missing_var_df = pd.DataFrame(missing_var_data)
        missing_var_df = missing_var_df.set_index("var_id")
        missing_var_df.index.name = None

        # Create zero X matrix for missing features
        x_dtype = adata.X.dtype if adata.X is not None else np.float32  # type: ignore[union-attr]
        if sp.issparse(adata.X):
            missing_x = sp.csr_matrix((n_cells, n_missing), dtype=x_dtype)
        else:
            missing_x = np.zeros((n_cells, n_missing), dtype=x_dtype)

        # Create AnnData for missing features
        missing_adata = AnnData(
            X=missing_x,
            obs=adata.obs.copy(),
            var=missing_var_df,
        )

        # Concatenate along var axis
        combined = concat([adata, missing_adata], axis=1, merge="first")
    else:
        combined = adata

    # Reorder to match schema feature order
    combined = combined[:, all_features].copy()  # type: ignore[index]

    # Update adata in-place
    adata._init_as_actual(
        X=combined.X,
        obs=combined.obs,
        var=combined.var,
        uns=combined.uns,
        obsm=combined.obsm,
        varm=combined.varm,
        obsp=combined.obsp,
        varp=combined.varp,
        layers=combined.layers,
        raw=combined.raw,
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

    Perform standard sanitization, then expand the var DataFrame to include
    all features from the ingest schema. Missing features are added with
    zero-filled X rows and var columns defaulted to 0. Features are reordered
    to match the schema feature order.

    :param adata: The AnnData object to sanitize. Modified in-place.
    :param ingest_schema: Schema containing the full feature set and var columns.
    """
    sanitize_for_ingest(adata=adata)
    _expand_to_full_feature_schema(adata=adata, ingest_schema=ingest_schema)
