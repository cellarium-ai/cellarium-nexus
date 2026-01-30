"""Kubeflow pipeline components for Nexus data operations."""

from cellarium.nexus.workflows.kubeflow.components.bq import (
    create_ingest_files_job,
    extract_job,
    ingest_data_to_bigquery_job,
    mark_curriculum_as_finished_job,
    prepare_extract_job,
)
from cellarium.nexus.workflows.kubeflow.components.soma import (
    mark_soma_curriculum_as_finished_job,
    soma_grouped_extract_job,
    soma_ingest_partition_job,
    soma_prepare_ingest_plan_job,
    soma_randomized_extract_job,
    soma_validate_sanitize_job,
)
from cellarium.nexus.workflows.kubeflow.components.validation import validate_anndata_files_job

__all__ = [
    # BQ components
    "create_ingest_files_job",
    "ingest_data_to_bigquery_job",
    "prepare_extract_job",
    "extract_job",
    "mark_curriculum_as_finished_job",
    # SOMA components
    "soma_randomized_extract_job",
    "soma_grouped_extract_job",
    "mark_soma_curriculum_as_finished_job",
    "soma_validate_sanitize_job",
    "soma_prepare_ingest_plan_job",
    "soma_ingest_partition_job",
    # Validation components
    "validate_anndata_files_job",
]
