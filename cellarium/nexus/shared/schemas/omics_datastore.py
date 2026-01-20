"""
Schema definitions for omics datastore.
"""

from __future__ import annotations

import base64
import io
from typing import Any, Literal

import pandas as pd
from pydantic import BaseModel, Field, PrivateAttr

# TileDB-SOMA supported dtype literals for obs/var columns
# Includes: primitive types, string, categorical, and timestamps
SomaDtype = Literal[
    # Boolean
    "bool",
    # Signed integers
    "int8",
    "int16",
    "int32",
    "int64",
    # Unsigned integers
    "uint8",
    "uint16",
    "uint32",
    "uint64",
    # Floating point
    "float32",
    "float64",
    # String types
    "str",
    "string",
    # Categorical
    "category",
    # Timestamps
    "datetime64[s]",
    "datetime64[ms]",
    "datetime64[us]",
    "datetime64[ns]",
]


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


class ObsSchemaDescriptor(BaseModel):
    """Descriptor for a single obs column validation."""

    name: str
    dtype: SomaDtype
    nullable: bool = False


class ExperimentVarSchema(BaseModel):
    """Complete var schema with canonical feature data.

    Store var DataFrame as base64-encoded Parquet for compactness.
    The DataFrame index contains feature IDs (e.g., ENSEMBL gene IDs).
    Used for first AnnData to define experiment var DataFrame, and for
    subsequent AnnDatas to validate feature subset membership.

    :param is_subset: If True, input AnnData may have fewer features than schema
        (e.g., scRNA-seq with GENCODE47 where file has 10k measured genes).
        If False, input must have exactly the same features (e.g., OPS imaging).
        Extra features not in schema are always rejected.
    """

    is_subset: bool = True
    _var_parquet_b64: str = PrivateAttr()
    _cached_df: pd.DataFrame | None = PrivateAttr(default=None)

    def __init__(self, *, _var_parquet_b64: str, is_subset: bool = True, **data: Any) -> None:
        """
        Initialize ExperimentVarSchema.

        :param _var_parquet_b64: Base64-encoded Parquet bytes of var DataFrame.
        :param is_subset: Whether input AnnData may have subset of schema features.
        """
        super().__init__(is_subset=is_subset, **data)
        self._var_parquet_b64 = _var_parquet_b64

    def to_dataframe(self) -> pd.DataFrame:
        """
        Deserialize var Parquet into DataFrame.

        :return: Var DataFrame with feature IDs as index.
        """
        if self._cached_df is None:
            parquet_bytes = base64.b64decode(self._var_parquet_b64)
            self._cached_df = pd.read_parquet(io.BytesIO(parquet_bytes))
        return self._cached_df.copy()

    def get_feature_ids(self) -> list[str]:
        """
        Extract ordered feature IDs from var index.

        :return: List of feature IDs in schema order.
        """
        return self.to_dataframe().index.tolist()

    @classmethod
    def from_dataframe(cls, *, var_df: pd.DataFrame, is_subset: bool = True) -> ExperimentVarSchema:
        """
        Create schema from pandas DataFrame.

        The DataFrame index must contain feature IDs (e.g., ENSEMBL gene IDs).
        All columns will be preserved as var metadata.

        :param var_df: Var DataFrame with feature IDs as index.
        :param is_subset: Whether input AnnData may have subset of schema features.

        :return: ExperimentVarSchema instance.
        """
        buffer = io.BytesIO()
        var_df.to_parquet(buffer, index=True)
        b64 = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return cls(_var_parquet_b64=b64, is_subset=is_subset)

    def model_dump(self, **kwargs: Any) -> dict[str, Any]:
        """
        Serialize to dict including private attribute.

        :return: Dictionary with _var_parquet_b64 field.
        """
        data = super().model_dump(**kwargs)
        data["_var_parquet_b64"] = self._var_parquet_b64
        return data

    @classmethod
    def model_validate(cls, obj: Any, **kwargs: Any) -> ExperimentVarSchema:
        """
        Deserialize from dict including private attribute.

        :param obj: Dictionary with _var_parquet_b64 field.

        :return: ExperimentVarSchema instance.
        """
        if isinstance(obj, dict) and "_var_parquet_b64" in obj:
            is_subset = obj.get("is_subset", True)
            return cls(_var_parquet_b64=obj["_var_parquet_b64"], is_subset=is_subset)
        return super().model_validate(obj, **kwargs)


class IngestSchema(BaseModel):
    """Schema for validating AnnData before SOMA ingest."""

    obs_columns: list[ObsSchemaDescriptor]
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
