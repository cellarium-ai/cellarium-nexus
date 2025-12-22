"""Kubeflow pipelines for Nexus data operations."""

from cellarium.nexus.workflows.kubeflow.pipelines.bq import (
    create_ingest_files_parallel_pipeline,
    extract_data_pipeline,
    ingest_data_to_bigquery_pipeline,
    prepare_extract_pipeline,
    run_extracts_pipeline,
)
from cellarium.nexus.workflows.kubeflow.pipelines.soma import (
    run_soma_grouped_extracts_pipeline,
    run_soma_randomized_extracts_pipeline,
    soma_grouped_extract_data_pipeline,
    soma_randomized_extract_data_pipeline,
)
from cellarium.nexus.workflows.kubeflow.pipelines.validation import validate_anndata_pipeline

__all__ = [
    # BQ pipelines
    "prepare_extract_pipeline",
    "run_extracts_pipeline",
    "create_ingest_files_parallel_pipeline",
    "ingest_data_to_bigquery_pipeline",
    "extract_data_pipeline",
    # SOMA pipelines
    "run_soma_randomized_extracts_pipeline",
    "run_soma_grouped_extracts_pipeline",
    "soma_randomized_extract_data_pipeline",
    "soma_grouped_extract_data_pipeline",
    # Validation pipelines
    "validate_anndata_pipeline",
]
