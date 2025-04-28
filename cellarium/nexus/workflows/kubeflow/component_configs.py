from typing import Any, Dict, List, Optional

from pydantic import BaseModel

from cellarium.nexus.shared import schemas


class IngestTaskConfig(BaseModel):
    """
    Combined configuration for a single create_ingest_files and
    ingest_data_to_bigquery task pair.

    :param project_id: GCP project ID
    :param nexus_backend_api_url: URL of the Nexus backend API
    :param bigquery_dataset: BigQuery dataset name
    :param data_source_path: Path to the input file for creation
    :param bucket_name: GCS bucket name for staging files
    :param ingest_bucket_path: Specific path within bucket for staging files for this task
    :param tag: Tag associated with the data source
    :param metadata_columns: Optional dictionary containing obs and var column mappings
    """

    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    data_source_path: str
    bucket_name: str
    ingest_bucket_path: str
    tag: str
    metadata_columns: Optional[Dict[str, Any]] = None


class BQOpsPrepareExtract(BaseModel):
    name: str
    creator_id: int
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    features: list[schemas.FeatureSchema]
    filters: dict[str, Any]
    obs_columns: list[str]
    extract_bin_size: int
    bucket_name: str
    extract_bucket_path: str


class ValidationConfig(BaseModel):
    """
    Configuration for validating AnnData files and reporting results.

    :param nexus_backend_api_url: URL of the Nexus backend API
    :param validation_report_id: ID of the validation report to update
    :param adata_gcs_paths: List of GCS paths to AnnData files to validate
    :param validation_methods: List of validation method names to apply to each file
    :param max_bytes_valid_per_file: Maximum file size in bytes for validation
    """

    nexus_backend_api_url: str
    validation_report_id: int
    adata_gcs_paths: list[str]
    validation_methods: list[str]
    max_bytes_valid_per_file: int


class BQOpsExtract(BaseModel):
    """
    Configuration for the extract job.

    :param project_id: GCP project ID
    :param nexus_backend_api_url: URL of the Nexus backend API
    :param bigquery_dataset: BigQuery dataset name
    :param name: Prefix for extract table names
    :param bins: List of bin numbers to extract
    :param bucket_name: GCS bucket name for output files
    :param extract_bucket_path: Path within bucket for output files
    :param obs_columns: Optional list of observation columns to include
    :param max_workers: Optional maximum number of parallel workers
    """

    name: str
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    bins: list[int]
    bucket_name: str
    extract_bucket_path: str
    obs_columns: list[str] | None = None
    max_workers: int | None = None
