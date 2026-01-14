"""
Sanitization functions for AnnData before SOMA ingest.

This module provides functions to remove unsupported AnnData slots
and prepare data for TileDB SOMA ingestion.
"""

from anndata import AnnData


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
