"""Coordinate AnnData preprocessing operations (validation and sanitization)."""

from __future__ import annotations

import logging

import anndata
from pydantic import ValidationError as PydanticValidationError

from cellarium.nexus.clients.nexus_backend_client import NexusBackendAPIClient
from cellarium.nexus.omics_datastore.soma_ops import TileDBSOMAIngestor
from cellarium.nexus.shared.schemas.omics_datastore import IngestSchema
from cellarium.nexus.shared.utils import WorkspaceFileManager, gcp

logger = logging.getLogger(__name__)


class PreprocessingCoordinator:
    """
    Coordinate AnnData preprocessing operations (validation and sanitization).

    This coordinator handles validation and sanitization workflows without
    requiring a SOMA experiment URI. For actual ingestion, use SomaIngestCoordinator.
    """

    def __init__(self, *, nexus_backend_api_url: str | None = None) -> None:
        """
        Initialize the preprocessing coordinator.

        :param nexus_backend_api_url: Optional URL for Nexus backend API (required for validation reporting)
        """
        self.ingestor = TileDBSOMAIngestor()
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url) if nexus_backend_api_url else None

        if nexus_backend_api_url is None:
            logger.info("Initialized PreprocessingCoordinator without backend API (reporting disabled)")

    @staticmethod
    def _parse_gcs_uri(gcs_uri: str) -> tuple[str, str]:
        """
        Parse GCS URI into bucket name and blob path.

        :param gcs_uri: GCS URI (e.g., gs://bucket/path/to/file)

        :return: Tuple of (bucket_name, blob_path)

        :raises ValueError: If GCS URI is invalid
        """
        bucket_name, blob_path = gcp.get_bucket_name_and_file_path_from_gc_path(full_gs_path=gcs_uri)
        if not bucket_name or not blob_path:
            raise ValueError(f"Invalid GCS URI: {gcs_uri}")
        return bucket_name, blob_path

    def validate_and_sanitize_files(
        self,
        *,
        ingest_schema_uri: str,
        input_h5ad_uris: list[str],
        output_h5ad_uris: list[str],
        validation_report_id: int | None = None,
    ) -> None:
        """
        Validate and sanitize AnnData files and optionally report results.

        Validate and sanitize each file, and create validation report items if backend API is available.
        Per-file semantics: is_valid=True only if validation AND sanitization both succeeded.

        :param ingest_schema_uri: GCS URI to the ingest schema JSON file
        :param input_h5ad_uris: List of GCS URIs for input h5ad files
        :param output_h5ad_uris: List of GCS URIs for sanitized output files
        :param validation_report_id: Optional ID of ValidationReport to attach items to (requires backend API)

        :raises IOError: If GCS operations fail
        :raises ValueError: If input/output URIs have mismatched lengths
        """
        if len(input_h5ad_uris) != len(output_h5ad_uris):
            raise ValueError("input_h5ad_uris and output_h5ad_uris must be the same length")

        if validation_report_id is not None and not self.backend_client:
            logger.warning("Validation report ID provided but backend API not configured - reporting disabled")

        logger.info(
            f"Starting validation (reporting {'enabled' if self.backend_client and validation_report_id else 'disabled'})"
        )

        schema_bucket, schema_blob = gcp.get_bucket_name_and_file_path_from_gc_path(ingest_schema_uri)
        schema_manager = WorkspaceFileManager(bucket_name=schema_bucket)
        ingest_schema_data = schema_manager.load_json_from_bucket(remote_path=schema_blob)
        ingest_schema = IngestSchema.model_validate(ingest_schema_data)

        file_results = []
        first_bucket, _ = self._parse_gcs_uri(input_h5ad_uris[0])
        workspace_manager = WorkspaceFileManager(bucket_name=first_bucket)

        with workspace_manager.temp_workspace() as workspace:
            for idx, (input_uri, output_uri) in enumerate(zip(input_h5ad_uris, output_h5ad_uris)):
                file_result = {
                    "input_uri": input_uri,
                    "output_uri": output_uri,
                    "is_valid": False,
                    "message": None,
                }

                logger.info(f"Validating and sanitizing {input_uri} -> {output_uri}")

                input_bucket, input_path = self._parse_gcs_uri(input_uri)
                output_bucket, output_path = self._parse_gcs_uri(output_uri)

                local_input = workspace["input"] / f"input_{idx}.h5ad"
                local_output = workspace["output"] / f"sanitized_{idx}.h5ad"

                input_manager = WorkspaceFileManager(bucket_name=input_bucket)
                input_manager.download_file_from_bucket(
                    remote_path=input_path,
                    local_path=local_input,
                )

                adata = anndata.read_h5ad(filename=local_input)

                try:
                    self.ingestor.validate_and_sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)
                except PydanticValidationError as e:
                    error_msg = f"Validation failed: {str(e)}"
                    file_result["message"] = error_msg
                    logger.warning(f"{input_uri}: {error_msg}")
                else:
                    adata.write_h5ad(filename=local_output)

                    output_manager = WorkspaceFileManager(bucket_name=output_bucket)
                    output_manager.upload_file_to_bucket(
                        local_path=local_output,
                        remote_path=output_path,
                    )

                    file_result["is_valid"] = True
                    logger.info(f"Successfully validated and sanitized {input_uri}")

                file_results.append(file_result)

                if self.backend_client and validation_report_id is not None:
                    sanitized_path = output_uri if file_result["is_valid"] else None
                    self.backend_client.create_validation_report_item(
                        report_id=validation_report_id,
                        input_file_path=input_uri,
                        validator_name="nexus.soma_ops.validate_and_sanitize",
                        is_valid=file_result["is_valid"],
                        message=file_result["message"],
                        sanitized_file_path=sanitized_path,
                    )
