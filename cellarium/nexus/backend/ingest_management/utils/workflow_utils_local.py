"""
Local workflow utilities for SOMA validation and sanitization pipeline.

Handles dataset schema conversion and validation report creation.
"""

from __future__ import annotations

import logging

from cellarium.nexus.backend.cell_management.models import OmicsDataset
from cellarium.nexus.backend.ingest_management.models import Ingest, IngestSchema
from cellarium.nexus.backend.ingest_management.utils.soma_schema_utils import django_ingest_schema_to_pydantic
from cellarium.nexus.coordinator.preprocessing_coordinator import PreprocessingCoordinator
from cellarium.nexus.coordinator.soma_ingest_coordinator import SomaIngestCoordinator

logger = logging.getLogger(__name__)


def run_validate_and_sanitize(
    *,
    ingest_schema: IngestSchema,
    input_h5ad_uris: list[str],
    output_h5ad_uris: list[str],
    validation_report_id: int | None = None,
    nexus_backend_api_url: str | None = None,
) -> None:
    """
    Validate and sanitize AnnData files against an ingest schema with optional reporting.

    Convert IngestSchema to Pydantic and run PreprocessingCoordinator.
    This is the main entry point for validation pipelines from the admin layer.

    :param ingest_schema: IngestSchema instance to validate against
    :param input_h5ad_uris: List of GCS URIs for input h5ad files
    :param output_h5ad_uris: List of GCS URIs for sanitized output files
    :param validation_report_id: Optional ID of ValidationReport to attach items to
    :param nexus_backend_api_url: Optional URL for Nexus backend API (required if validation_report_id provided)

    :raises ValueError: If schema is invalid
    :raises IOError: If GCS operations fail
    """
    logger.info(f"Loading schema '{ingest_schema.name}' for validation")
    pydantic_schema = django_ingest_schema_to_pydantic(django_schema=ingest_schema)

    coordinator = PreprocessingCoordinator(nexus_backend_api_url=nexus_backend_api_url)

    logger.info(
        f"Starting validation for schema '{ingest_schema.name}' "
        f"(reporting {'enabled' if validation_report_id else 'disabled'})"
    )

    coordinator.validate_and_sanitize_files(
        ingest_schema=pydantic_schema,
        input_h5ad_uris=input_h5ad_uris,
        output_h5ad_uris=output_h5ad_uris,
        validation_report_id=validation_report_id,
    )


def run_soma_ingest(
    *,
    omics_dataset: OmicsDataset,
    ingest: Ingest,
    sanitized_h5ad_uris: list[str],
    ingest_batch_size: int,
    measurement_name: str,
    nexus_backend_api_url: str | None = None,
) -> None:
    """
    Ingest sanitized h5ad files into a SOMA experiment.

    Load dataset schema, prepare ingest plan, and sequentially ingest partitions.
    This is the main entry point for ingest operations from the admin layer.

    :param omics_dataset: OmicsDataset instance to ingest into
    :param ingest: Ingest instance for tracking progress
    :param sanitized_h5ad_uris: List of GCS URIs for sanitized h5ad files
    :param ingest_batch_size: Number of h5ad files per partition
    :param measurement_name: Name of the SOMA measurement (e.g., RNA, ATAC)
    :param nexus_backend_api_url: Optional URL for Nexus backend API

    :raises ValueError: If dataset not found or has no schema
    :raises IOError: If GCS operations fail
    """
    if not sanitized_h5ad_uris:
        raise ValueError("No h5ad files provided for ingestion")

    if omics_dataset.schema is None:
        raise ValueError(f"Dataset '{omics_dataset.name}' has no attached schema. Ingest requires a schema.")

    pydantic_schema = django_ingest_schema_to_pydantic(django_schema=omics_dataset.schema)

    coordinator = SomaIngestCoordinator(experiment_uri=omics_dataset.uri, nexus_backend_api_url=nexus_backend_api_url)

    # Prepare ingest plan
    ingest_plan = coordinator.prepare_ingest_plan(
        h5ad_uris=sanitized_h5ad_uris,
        measurement_name=measurement_name,
        ingest_schema=pydantic_schema,
        ingest_batch_size=ingest_batch_size,
        first_adata_gcs_path=sanitized_h5ad_uris[0],
    )

    # Ingest each partition sequentially
    num_partitions = (ingest_plan.total_files + ingest_batch_size - 1) // ingest_batch_size
    for partition_index in range(num_partitions):
        logger.info(f"Ingesting partition {partition_index + 1}/{num_partitions}")
        coordinator.ingest_partition(
            ingest_plan=ingest_plan,
            partition_index=partition_index,
            h5ad_file_paths=sanitized_h5ad_uris,
            omics_dataset_name=omics_dataset.name,
            parent_ingest_id=ingest.id,
        )
