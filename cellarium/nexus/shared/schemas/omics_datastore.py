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
    :param measured_genes_mask: Matrix indicating which genes are measured for each tag

    :raise ValueError: If validation fails
    """

    total_bins: int
    last_bin_size: int
    total_cells: int = 0
    filters: dict[str, Any] | None = None
    extract_bin_size: int | None = None
    category_metadata: dict[str, list[str]] = Field(default_factory=dict)
    measured_genes_mask: list[dict[str, Any]] = Field(default_factory=list)

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
