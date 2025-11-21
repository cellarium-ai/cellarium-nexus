"""
TileDB SOMA operations for the Nexus omics datastore.

This package provides extract operations for TileDB SOMA experiments.
"""

from cellarium.nexus.omics_datastore.soma_ops.data_operator import TileDBSOMADataOperator
from cellarium.nexus.omics_datastore.soma_ops.exceptions import (
    SomaExtractError,
    SomaFilterError,
    SomaOperationError,
    SomaPlanningError,
    SomaReadError,
    SomaWriteError,
)
from cellarium.nexus.omics_datastore.soma_ops.extract import shuffle_extracted_chunks
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter
from cellarium.nexus.omics_datastore.soma_ops.planning import plan_soma_extract

__all__ = [
    "TileDBSOMADataOperator",
    "build_soma_value_filter",
    "plan_soma_extract",
    "shuffle_extracted_chunks",
    # Exceptions
    "SomaOperationError",
    "SomaReadError",
    "SomaWriteError",
    "SomaFilterError",
    "SomaExtractError",
    "SomaPlanningError",
]
