"""
TileDB SOMA operations for the Nexus omics datastore.

This package provides extract operations for TileDB SOMA experiments.
"""

from cellarium.nexus.omics_datastore.soma_ops.data_ingestor import TileDBSOMAIngestor
from cellarium.nexus.omics_datastore.soma_ops.data_operator import TileDBSOMADataOperator
from cellarium.nexus.omics_datastore.soma_ops.exceptions import (
    SomaExtractError,
    SomaFilterError,
    SomaOperationError,
    SomaPrepareCurriculumMetadataError,
    SomaReadError,
    SomaWriteError,
)

__all__ = [
    "TileDBSOMAIngestor",
    "TileDBSOMADataOperator",
    # Exceptions
    "SomaOperationError",
    "SomaReadError",
    "SomaWriteError",
    "SomaFilterError",
    "SomaExtractError",
    "SomaPrepareCurriculumMetadataError",
]
