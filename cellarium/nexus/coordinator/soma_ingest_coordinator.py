"""
Coordinate TileDB SOMA ingest workflows for Nexus.
"""

from __future__ import annotations

import logging

import anndata
from pydantic import ValidationError as PydanticValidationError

from cellarium.nexus.clients.nexus_backend_client import NexusBackendAPIClient
from cellarium.nexus.omics_datastore.soma_ops import TileDBSOMAIngestor
from cellarium.nexus.omics_datastore.soma_ops.utils import get_block_slice
from cellarium.nexus.shared.schemas.omics_datastore import IngestPlanMetadata, IngestSchema
from cellarium.nexus.shared.utils import WorkspaceFileManager, gcp

logger = logging.getLogger(__name__)


class SomaIngestCoordinator:
    """
    Coordinate TileDB SOMA ingest operations.

    Provides high-level helpers for:
    - validating and sanitizing input AnnData files
    - preparing ingest plans
    - ingesting a partition of files into a SOMA experiment
    """

    def __init__(self, *, experiment_uri: str, nexus_backend_api_url: str | None = None) -> None:
        """
        Initialize the SOMA ingest coordinator.

        :param experiment_uri: URI of the SOMA experiment
        :param nexus_backend_api_url: Optional URL for Nexus backend API (required for validation reporting)

        :raise ValueError: If experiment_uri is empty
        """
        if not experiment_uri:
            raise ValueError("experiment_uri cannot be empty")

        self.experiment_uri = experiment_uri
        self.ingestor = TileDBSOMAIngestor()
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url) if nexus_backend_api_url else None

        if nexus_backend_api_url is None:
            logger.info("Initialized SomaIngestCoordinator without backend API (reporting disabled)")

        logger.info(f"Experiment URI {experiment_uri}")

    @staticmethod
    def _parse_gcs_uri(gcs_uri: str) -> tuple[str, str]:
        bucket_name, blob_path = gcp.get_bucket_name_and_file_path_from_gc_path(full_gs_path=gcs_uri)
        if not bucket_name or not blob_path:
            raise ValueError(f"Invalid GCS URI: {gcs_uri}")
        return bucket_name, blob_path

    def validate_and_sanitize_files(
        self,
        *,
        ingest_schema: IngestSchema,
        input_h5ad_uris: list[str],
        output_h5ad_uris: list[str],
        validation_report_id: int | None = None,
    ) -> None:
        """
        Validate and sanitize AnnData files and optionally report results.

        Validate and sanitize each file, and create validation report items if backend API is available.
        Per-file semantics: is_valid=True only if validation AND sanitization both succeeded.

        :param ingest_schema: Pydantic schema for validation and sanitization
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
                    self.backend_client.create_validation_report_item(
                        report_id=validation_report_id,
                        input_file_gcs_path=input_uri,
                        validator_name="nexus.soma_ops.validate_and_sanitize",
                        is_valid=file_result["is_valid"],
                        message=file_result["message"],
                    )

    def prepare_ingest_plan(
        self,
        *,
        h5ad_uris: list[str],
        measurement_name: str,
        ingest_schema: IngestSchema,
        ingest_batch_size: int,
        first_adata_gcs_path: str | None = None,
    ) -> IngestPlanMetadata:
        """
        Prepare a partitioned SOMA ingest plan.

        :param h5ad_uris: List of GCS URIs for h5ad files to ingest
        :param measurement_name: Name of the SOMA measurement
        :param ingest_schema: Schema for validation and sanitization
        :param ingest_batch_size: Number of h5ads per partition
        :param first_adata_gcs_path: Optional GCS URI for first AnnData used to create schema-only experiment

        :raise ValueError: If ingest plan preparation fails

        :return: IngestPlanMetadata
        """
        first_adata = None
        if first_adata_gcs_path:
            logger.info(f"Loading first AnnData from {first_adata_gcs_path}")
            bucket_name, blob_path = self._parse_gcs_uri(first_adata_gcs_path)
            workspace_manager = WorkspaceFileManager(bucket_name=bucket_name)
            with workspace_manager.temp_workspace() as workspace:
                local_path = workspace["input"] / "first_adata.h5ad"
                workspace_manager.download_file_from_bucket(
                    remote_path=blob_path,
                    local_path=local_path,
                )
                # Take the first 10 cells to avoid loading the entire dataset into memory.
                # The ingestor only needs anndata for using it as a schema to create the experiment
                first_adata = anndata.read_h5ad(filename=local_path)[:10]

        return self.ingestor.prepare_ingest_plan(
            experiment_uri=self.experiment_uri,
            h5ad_paths=h5ad_uris,
            measurement_name=measurement_name,
            ingest_schema=ingest_schema,
            ingest_batch_size=ingest_batch_size,
            first_adata=first_adata,
        )

    def save_ingest_plan_to_gcs(self, *, ingest_plan: IngestPlanMetadata, ingest_plan_gcs_path: str) -> str:
        """
        Serialize and save an ingest plan to GCS.

        :param ingest_plan: Ingest plan metadata
        :param ingest_plan_gcs_path: GCS URI where plan should be stored

        :return: Full GCS URI to the stored plan
        """
        bucket_name, blob_path = self._parse_gcs_uri(ingest_plan_gcs_path)
        workspace_manager = WorkspaceFileManager(bucket_name=bucket_name)
        return workspace_manager.save_json_to_bucket(
            data=ingest_plan.model_dump(),
            remote_path=blob_path,
        )

    def load_ingest_plan_from_gcs(self, *, ingest_plan_gcs_path: str) -> IngestPlanMetadata:
        """
        Load and deserialize an ingest plan from GCS.

        :param ingest_plan_gcs_path: GCS URI for the ingest plan JSON

        :return: IngestPlanMetadata
        """
        bucket_name, blob_path = self._parse_gcs_uri(ingest_plan_gcs_path)
        workspace_manager = WorkspaceFileManager(bucket_name=bucket_name)
        data = workspace_manager.load_json_from_bucket(remote_path=blob_path)
        return IngestPlanMetadata.model_validate(obj=data)

    def ingest_partition(
        self,
        *,
        ingest_plan: IngestPlanMetadata,
        partition_index: int,
        h5ad_uris: list[str] | None = None,
    ) -> None:
        """
        Ingest a single partition of h5ad files into a SOMA experiment.

        :param ingest_plan: Ingest plan metadata
        :param partition_index: Zero-based partition index
        :param h5ad_uris: Optional list of GCS URIs to use instead of plan sources

        :raise ValueError: If the provided URIs do not match expected partition size
        """
        source_uris = h5ad_uris if h5ad_uris is not None else ingest_plan.source_h5ad_uris
        slice_start, slice_end = get_block_slice(
            total_items=ingest_plan.total_files,
            partition_index=partition_index,
            block_size=ingest_plan.ingest_batch_size,
        )
        partition_uris = source_uris[slice_start:slice_end]

        if not partition_uris:
            logger.info(f"Partition {partition_index} has no files to ingest")
            return

        first_bucket, _ = self._parse_gcs_uri(partition_uris[0])
        workspace_manager = WorkspaceFileManager(bucket_name=first_bucket)
        with workspace_manager.temp_workspace() as workspace:
            local_paths: list[str] = []

            for idx, uri in enumerate(partition_uris):
                bucket_name, blob_path = self._parse_gcs_uri(uri)
                local_path = workspace["input"] / f"ingest_{partition_index}_{idx}.h5ad"
                bucket_manager = WorkspaceFileManager(bucket_name=bucket_name)
                bucket_manager.download_file_from_bucket(
                    remote_path=blob_path,
                    local_path=local_path,
                )
                local_paths.append(str(local_path))

            self.ingestor.ingest_h5ads_partition(
                ingest_plan=ingest_plan,
                partition_index=partition_index,
                local_h5ad_paths=local_paths,
            )
