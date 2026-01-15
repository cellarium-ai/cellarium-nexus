"""
Schema definitions for omics datastore.
"""

from typing import Any, Literal

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


class GroupedBin(BaseModel):
    """A bin containing cells from a specific group combination."""

    group_key: str
    group_filter: str
    joinid_min: int
    joinid_max: int
    cell_count: int


class BaseCurriculumMetadata(BaseModel):
    """Base class for SOMA curriculum metadata."""

    experiment_uri: str
    value_filter: str
    total_cells: int
    num_bins: int
    extract_bin_size: int
    last_bin_size: int
    filters: dict[str, Any] | None = None
    var_joinids: list[int] | None = None
    var_filter_column: str | None = None
    var_filter_values: list[str] | None = None
    obs_columns: list[str] | None = None
    var_columns: list[str] | None = None
    x_layer: str = "X"


class RandomizedCurriculumMetadata(BaseCurriculumMetadata):
    """Curriculum metadata for randomized extraction with contiguous ranges."""

    id_ranges: list[IdContiguousRange]
    range_size: int
    num_ranges: int
    extract_bin_indexes: list[int]


class GroupedCurriculumMetadata(BaseCurriculumMetadata):
    """Curriculum metadata for grouped extraction."""

    extract_bin_keys: list[str]
    grouped_bins: list[GroupedBin]


# Ingest validation schemas


class ObsDescriptor(BaseModel):
    """Descriptor for a single obs column validation."""

    name: str
    dtype: str
    nullable: bool = False


class VarDescriptor(BaseModel):
    """Descriptor for a single var column."""

    name: str
    dtype: str
    nullable: bool = False


class ExperimentVarSchema(BaseModel):
    """Schema for var (features) validation.

    Define only feature ids without additional descriptors.
    """

    features: list[str]
    is_subset: bool = True


class IngestSchema(BaseModel):
    """Schema for validating AnnData before SOMA ingest."""

    obs_columns: list[ObsDescriptor]
    var_schema: ExperimentVarSchema
    x_validation_type: Literal["count_matrix", "feature_matrix"]


class IngestPlanMetadata(BaseModel):
    """
    Metadata for a partitioned SOMA ingest plan.

    This schema captures all information needed to execute a partitioned ingest
    workflow. The `source_h5ad_uris` field contains remote URIs (e.g., GCS paths)
    representing the logical dataset. The actual ingest function receives local
    file paths after the coordinator downloads the files.

    The `registration_mapping` field contains a base64-encoded pickle of
    `tiledbsoma.io.ExperimentAmbientLabelMapping` that can be deserialized
    for use during ingestion.
    """

    experiment_uri: str
    source_h5ad_uris: list[str]
    measurement_name: str
    total_files: int
    ingest_batch_size: int
    num_partitions: int
    last_partition_size: int
    ingest_schema: IngestSchema
    registration_mapping_pickle: str
