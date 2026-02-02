from datetime import datetime
from typing import Any

from pydantic import BaseModel


class IngestInfoAPISchema(BaseModel):
    id: int
    nexus_uuid: str
    omics_dataset: str | None = None
    ingest_id: int | None = None
    ingest_start_timestamp: datetime
    ingest_finish_timestamp: datetime | None = None
    gcs_file_path: str | None = None
    tag: str | None = None
    metadata_extra: dict[str, Any] | None = None
    status: str | None = None


class FeatureInfoAPISchema(BaseModel):
    id: int
    ensemble_id: str
    symbol: str
    ingest_id: int
    tag: str | None = None
    metadata_extra: dict[str, Any] | None = None


class CurriculumAPISchema(BaseModel):
    id: int
    creator_id: int
    created_at: datetime
    cell_count: int | None = None
    name: str | None = None
    description: str | None = None
    extract_files_path: str | None = None
    metadata_file_path: str | None = None
    extract_bin_size: int | None = None
    filters_json: dict[str, Any] | None = None


class ValidationReportAPISchema(BaseModel):
    id: int
    creator_id: int | None = None
    created_at: datetime


class ValidationReportItemAPISchema(BaseModel):
    id: int
    report_id: int
    input_file_gcs_path: str
    validator_name: str
    is_valid: bool
    message: str | None = None
    created_at: datetime


class BackendResetAPISchema(BaseModel):
    status: str
    message: str
    repopulated_keys: list[str]
    count: int
