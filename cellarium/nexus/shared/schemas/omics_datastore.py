"""
Schema definitions for omics datastore.
"""

from typing import Any

from pydantic import BaseModel, Field


class FeatureSchema(BaseModel):
    """
    Schema for feature data.

    :param id: Unique identifier for the feature
    :param symbol: Gene symbol
    :param ensemble_id: Ensemble identifier
    """

    id: int
    symbol: str
    ensemble_id: str


class ExtractMetadata(BaseModel):
    """
    Store metadata about the extract.

    :param total_bins: Total number of extract bins
    :param last_bin_size: Size of the last bin
    :param total_cells: Total number of cells in the extract
    :param filters: Filters used for the extract
    :param extract_bin_size: Size of extract bins
    :param category_metadata: Metadata about categorical columns

    :raise ValueError: If validation fails
    """

    total_bins: int
    last_bin_size: int
    total_cells: int = 0
    filters: dict[str, Any] | None = None
    extract_bin_size: int | None = None
    category_metadata: dict[str, list[str]] = Field(default_factory=dict)

    def calculate_cell_count(self) -> int:
        """
        Get the total cell count.

        If total_cells is set, returns that value. Otherwise, calculates based on bin sizes.

        :return: Total cell count
        """
        if self.total_cells > 0:
            return self.total_cells

        bin_size = self.extract_bin_size or 0
        return bin_size * (self.total_bins - 1) + self.last_bin_size if self.total_bins > 0 else 0


class SomaJoinIdRange(BaseModel):
    """
    Store a contiguous joinid range.

    :param start: Inclusive start soma_joinid
    :param end: Inclusive end soma_joinid
    """

    start: int
    end: int


class SomaExtractPlan(BaseModel):
    """
    Store SOMA extract plan information.

    :param experiment_uri: URI of the SOMA experiment
    :param value_filter: SOMA obs filter expression used for this plan
    :param joinid_ranges: List of contiguous soma_joinid ranges
    :param total_cells: Total number of cells matching the filter
    :param range_size: Target number of cells per range (for extraction)
    :param output_chunk_size: Target number of cells per output chunk (for shuffling)
    :param filters: Structured filter specification, using Nexus format
    :param var_joinids: Ordered list of feature soma_joinids to retain (optional)
    """

    experiment_uri: str
    value_filter: str
    joinid_ranges: list[SomaJoinIdRange]
    total_cells: int
    range_size: int
    output_chunk_size: int
    filters: dict[str, Any] | None = None
    var_joinids: list[int] | None = None
