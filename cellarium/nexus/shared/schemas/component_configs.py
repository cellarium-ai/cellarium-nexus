from typing import Any, Dict, Optional

from pydantic import BaseModel

from cellarium.nexus.shared.schemas.omics_datastore import FeatureSchema


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
