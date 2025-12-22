"""
TileDB SOMA operations for the Nexus omics datastore.

This package provides extract operations for TileDB SOMA experiments.
"""

from cellarium.nexus.omics_datastore.soma_ops.curriculum_grouped import prepare_grouped_curriculum
from cellarium.nexus.omics_datastore.soma_ops.curriculum_randomized import prepare_extract_curriculum
from cellarium.nexus.omics_datastore.soma_ops.data_operator import TileDBSOMADataOperator
from cellarium.nexus.omics_datastore.soma_ops.exceptions import (
    SomaExtractError,
    SomaFilterError,
    SomaOperationError,
    SomaPrepareCurriculumMetadataError,
    SomaReadError,
    SomaWriteError,
)
from cellarium.nexus.omics_datastore.soma_ops.extract_grouped import extract_grouped_bins
from cellarium.nexus.omics_datastore.soma_ops.extract_randomized import shuffle_extracted_chunks
from cellarium.nexus.omics_datastore.soma_ops.filters import build_soma_value_filter

__all__ = [
    "TileDBSOMADataOperator",
    "build_soma_value_filter",
    "prepare_extract_curriculum",
    "prepare_grouped_curriculum",
    "shuffle_extracted_chunks",
    "extract_grouped_bins",
    # Exceptions
    "SomaOperationError",
    "SomaReadError",
    "SomaWriteError",
    "SomaFilterError",
    "SomaExtractError",
    "SomaPrepareCurriculumMetadataError",
]
