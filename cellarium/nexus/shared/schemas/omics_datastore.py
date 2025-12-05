"""
Schema definitions for omics datastore.
"""

from typing import Any

from pydantic import BaseModel, Field


class FeatureSchema(BaseModel):
    id: int
    symbol: str
    ensemble_id: str


class ExtractMetadata(BaseModel):
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


class IdContiguousRange(BaseModel):
    start: int
    end: int


class SomaExtractPlan(BaseModel):
    experiment_uri: str
    value_filter: str
    id_ranges: list[IdContiguousRange]
    total_cells: int
    range_size: int
    num_output_chunks: int
    output_chunk_size: int
    output_chunk_indexes: list[int] | None = None
    filters: dict[str, Any] | None = None
    var_joinids: list[int] | None = None
    var_filter_column: str | None = None
    var_filter_values: list[str] | None = None
    obs_columns: list[str] | None = None
    var_columns: list[str] | None = None
    x_layer: str = "X"
