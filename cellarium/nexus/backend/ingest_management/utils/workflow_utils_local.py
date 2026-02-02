"""
Local workflow utilities for SOMA validation and sanitization pipeline.

Handles dataset schema conversion and validation report creation.
"""

from __future__ import annotations

import logging

from cellarium.nexus.backend.cell_management.models import OmicsDataset
from cellarium.nexus.backend.ingest_management.utils.soma_schema_utils import django_ingest_schema_to_pydantic
from cellarium.nexus.coordinator.soma_ingest_coordinator import SomaIngestCoordinator

logger = logging.getLogger(__name__)


def run_validate_and_sanitize(
    *,
    dataset_name: str,
    input_h5ad_uris: list[str],
    output_h5ad_uris: list[str],
    validation_report_id: int | None = None,
    experiment_uri: str | None = None,
    nexus_backend_api_url: str | None = None,
) -> None:
    """
    Validate and sanitize AnnData files for a dataset with optional reporting.

    Load dataset schema, convert to Pydantic, and run coordinator.
    This is the main entry point for validation pipelines from the admin layer.

    :param dataset_name: Name of the OmicsDataset to validate against
    :param input_h5ad_uris: List of GCS URIs for input h5ad files
    :param output_h5ad_uris: List of GCS URIs for sanitized output files
    :param validation_report_id: Optional ID of ValidationReport to attach items to
    :param experiment_uri: Optional URI of the target SOMA experiment (defaults to dataset experiment)
    :param nexus_backend_api_url: Optional URL for Nexus backend API (required if validation_report_id provided)

    :raises ValueError: If dataset not found or has no schema
    :raises IOError: If GCS operations fail
    """
    try:
        dataset = OmicsDataset.objects.get(name=dataset_name)
    except OmicsDataset.DoesNotExist:
        raise ValueError(f"Dataset '{dataset_name}' not found")

    if dataset.schema is None:
        raise ValueError(
            f"Dataset '{dataset_name}' has no attached schema. " "Validation and sanitization require a schema."
        )

    logger.info(f"Loading schema for dataset '{dataset_name}'")
    pydantic_schema = django_ingest_schema_to_pydantic(django_schema=dataset.schema)

    # Use dataset's experiment URI if not provided
    exp_uri = experiment_uri or dataset.uri
    if not exp_uri:
        raise ValueError(f"No experiment URI found for dataset '{dataset_name}'")

    coordinator = SomaIngestCoordinator(
        experiment_uri=exp_uri,
        nexus_backend_api_url=nexus_backend_api_url,
    )

    logger.info(
        f"Starting validation for dataset '{dataset_name}' "
        f"(reporting {'enabled' if validation_report_id else 'disabled'})"
    )

    coordinator.validate_and_sanitize_files(
        ingest_schema=pydantic_schema,
        input_h5ad_uris=input_h5ad_uris,
        output_h5ad_uris=output_h5ad_uris,
        validation_report_id=validation_report_id,
    )
