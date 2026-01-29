from typing import Any, Literal

from pydantic import BaseModel

from cellarium.nexus.shared.schemas.omics_datastore import FeatureSchema, IngestSchema


class CreateIngestFilesConfig(BaseModel):
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    input_file_path: str
    bucket_name: str
    bucket_stage_dir: str
    tag: str
    max_input_data_size: int
    column_mapping: dict[str, Any] | None = None
    validation_methods: list[str] | None = None
    uns_keys_to_keep: list[str] | None = None


class IngestFilesConfig(BaseModel):
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    bucket_name: str
    bucket_stage_dirs: list[str]
    num_workers: int = 4


class BQOpsPrepareExtract(BaseModel):
    extract_name: str
    creator_id: int
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    features: list[FeatureSchema]
    categorical_column_count_limit: int
    filters: dict[str, Any]
    obs_columns: list[str]
    extract_bin_size: int
    bucket_name: str
    extract_bucket_path: str
    extract_bin_keys: list[str] | None = None
    metadata_extra_columns: list[str] | None = None


class ValidationConfig(BaseModel):
    nexus_backend_api_url: str
    validation_report_id: int
    adata_gcs_paths: list[str]
    validation_methods: list[str]
    max_bytes_valid_per_file: int


class BQOpsExtract(BaseModel):
    extract_name: str
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    bins: list[int]
    bucket_name: str
    extract_bucket_path: str
    obs_columns: list[str] | None = None
    max_workers: int | None = None


class SomaOpsExtract(BaseModel):
    """Configuration for running a SOMA extract."""

    extract_name: str
    experiment_uri: str
    nexus_backend_api_url: str
    bucket_name: str
    extract_metadata_path: str
    extract_bucket_path: str
    partition_index: int = 0
    output_format: str = "h5ad"
    max_workers_extract: int | None = None

    # Discriminator field
    extract_type: Literal["randomized", "grouped"] = "randomized"

    # Randomized-specific params
    max_ranges_per_partition: int | None = None
    max_workers_shuffle: int | None = None

    # Grouped-specific params
    max_bins_per_partition: int | None = None


class SomaValidateSanitizeConfig(BaseModel):
    """Configuration for validating and sanitizing SOMA ingest inputs."""

    experiment_uri: str
    input_h5ad_uris: list[str]
    output_h5ad_uris: list[str]
    ingest_schema: IngestSchema
    max_bytes_per_file: int | None = None


class SomaIngestPlanConfig(BaseModel):
    """Configuration for preparing a SOMA ingest plan."""

    experiment_uri: str
    measurement_name: str
    ingest_schema: IngestSchema
    ingest_batch_size: int
    h5ad_uris: list[str]
    ingest_plan_gcs_path: str
    first_adata_gcs_path: str | None = None


class SomaIngestPartitionConfig(BaseModel):
    """Configuration for ingesting a SOMA partition."""

    experiment_uri: str
    ingest_plan_gcs_path: str
    partition_index: int
