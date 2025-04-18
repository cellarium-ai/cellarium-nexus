from pydantic import BaseModel
from typing import Dict, Any, Optional
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
    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    extract_table_prefix: str
    features: list[schemas.FeatureSchema]
    filters: dict[str, Any]
    obs_columns: list[str]
    extract_bin_size: int
    bucket_name: str
    extract_bucket_path: str


class BQOpsExtract(BaseModel):
    """
    Configuration for the extract job.

    :param project_id: GCP project ID
    :param nexus_backend_api_url: URL of the Nexus backend API
    :param bigquery_dataset: BigQuery dataset name
    :param extract_table_prefix: Prefix for extract table names
    :param bins: List of bin numbers to extract
    :param bucket_name: GCS bucket name for output files
    :param extract_bucket_path: Path within bucket for output files
    :param obs_columns: Optional list of observation columns to include
    :param max_workers: Optional maximum number of parallel workers
    """

    project_id: str
    nexus_backend_api_url: str
    bigquery_dataset: str
    extract_table_prefix: str
    bins: list[int]
    bucket_name: str
    extract_bucket_path: str
    obs_columns: list[str] | None = None
    max_workers: int | None = None
