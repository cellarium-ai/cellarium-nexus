"""
Coordinate TileDB SOMA ingest workflows for Nexus.
"""

from __future__ import annotations

import datetime
import logging

import anndata

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

    def _update_ingest_info_with_error(self, *, ingest_id: int, error_message: str | dict[str, str]) -> None:
        """
        Update ingest info with error message and mark as failed.

        :param ingest_id: ID of the ingest to update
        :param error_message: Error message to store
        """
        if not self.backend_client:
            logger.warning(f"Cannot report error for ingest {ingest_id}: backend client not configured")
            return

        try:
            self.backend_client.update_ingest_metadata_extra(
                ingest_id=ingest_id, new_metadata_extra={"error_occurred": error_message}
            )
            self.backend_client.update_ingest_status(ingest_id=ingest_id, new_status="FAILED")
        except Exception as e:
            logger.error(f"Failed to report error status for ingest {ingest_id}: {e}")

    def _update_ingest_info_with_success(self, *, ingest_id: int) -> None:
        """
        Mark ingest info as successfully completed.

        :param ingest_id: ID of the ingest to update
        """
        if not self.backend_client:
            logger.warning(f"Cannot report success for ingest {ingest_id}: backend client not configured")
            return

        try:
            self.backend_client.update_ingest_status(
                ingest_id=ingest_id, new_status="SUCCEEDED", ingest_finish_timestamp=datetime.datetime.now()
            )
        except Exception as e:
            logger.error(f"Failed to report success status for ingest {ingest_id}: {e}")

    def _create_ingest_info_records(
        self,
        *,
        partition_index: int,
        partition_uris: list[str],
        omics_dataset_name: str | None,
        parent_ingest_id: int | None,
    ) -> list[int]:
        """
        Create IngestInfo records for each file in the partition.

        :param partition_index: Zero-based partition index
        :param partition_uris: List of GCS URIs for this partition
        :param omics_dataset_name: Dataset name for creating records
        :param parent_ingest_id: Optional parent ingest ID

        :return: List of created IngestInfo IDs
        """
        ingest_file_ids: list[int] = []

        if not self.backend_client or not omics_dataset_name:
            if not omics_dataset_name:
                logger.info(
                    f"Partition {partition_index} ingest without backend reporting (omics_dataset_name not provided)"
                )
            return ingest_file_ids

        try:
            for uri in partition_uris:
                ingest_info = self.backend_client.create_ingest_file_info(
                    omics_dataset=omics_dataset_name, ingest_id=parent_ingest_id, gcs_file_path=uri
                )
                ingest_file_ids.append(ingest_info.id)
            logger.info(f"Created {len(ingest_file_ids)} IngestInfo records for partition {partition_index}")
        except Exception as e:
            logger.error(f"Failed to create IngestInfo records for partition {partition_index}: {e}")

        return ingest_file_ids

    def _report_ingest_started(self, *, ingest_file_ids: list[int]) -> None:
        """
        Report STARTED status for created IngestInfo records.

        :param ingest_file_ids: List of IngestInfo IDs to report
        """
        if not self.backend_client:
            return

        for ingest_id in ingest_file_ids:
            try:
                self.backend_client.update_ingest_status(ingest_id=ingest_id, new_status="STARTED")
            except Exception as e:
                logger.error(f"Failed to report STARTED status for ingest {ingest_id}: {e}")

    def _report_ingest_succeeded(self, *, ingest_file_ids: list[int], partition_index: int) -> None:
        """
        Report SUCCEEDED status for IngestInfo records.

        :param ingest_file_ids: List of IngestInfo IDs to report
        :param partition_index: Partition index for logging
        """
        for ingest_id in ingest_file_ids:
            self._update_ingest_info_with_success(ingest_id=ingest_id)
        logger.info(f"Partition {partition_index} ingested successfully, reported to backend")

    def _report_ingest_failed(
        self,
        *,
        ingest_file_ids: list[int],
        error_message: str,
        partition_index: int,
    ) -> None:
        """
        Report FAILED status for IngestInfo records.

        :param ingest_file_ids: List of IngestInfo IDs to report
        :param error_message: Error message to include
        :param partition_index: Partition index for logging
        """
        for ingest_id in ingest_file_ids:
            self._update_ingest_info_with_error(ingest_id=ingest_id, error_message=error_message)
        logger.error(f"Partition {partition_index} ingest failed, reported to backend: {error_message}")

    def _download_partition_files(
        self,
        *,
        partition_index: int,
        partition_uris: list[str],
        workspace: dict,
    ) -> list[str]:
        """
        Download h5ad files for a partition to local workspace.

        :param partition_index: Zero-based partition index
        :param partition_uris: List of GCS URIs for this partition
        :param workspace: Workspace dictionary with 'input' path

        :return: List of local file paths
        """
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

        return local_paths

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
        h5ad_file_paths: list[str],
        omics_dataset_name: str | None = None,
        parent_ingest_id: int | None = None,
    ) -> None:
        """
        Ingest a single partition of h5ad files into a SOMA experiment.

        :param ingest_plan: Ingest plan metadata
        :param partition_index: Zero-based partition index
        :param h5ad_file_paths: List of GCS URIs to ingest
        :param omics_dataset_name: Optional dataset name for creating IngestInfo records
        :param parent_ingest_id: Optional parent ingest ID

        :raise ValueError: If the provided URIs do not match expected partition size
        """
        dataset_name_or_ingest_id_missing = any([not omics_dataset_name, parent_ingest_id is None])
        if dataset_name_or_ingest_id_missing and self.backend_client is not None:
            raise ValueError(
                "Ingest plan metadata and parent ingest ID must be provided to report ingest status to backend. "
                "Either provide both or disable backend reporting by reinitializing the coordinator with "
                "backend_api_url=None."
            )
        slice_start, slice_end = get_block_slice(
            total_items=ingest_plan.total_files,
            partition_index=partition_index,
            block_size=ingest_plan.ingest_batch_size,
        )
        partition_uris = h5ad_file_paths[slice_start:slice_end]

        if not partition_uris:
            logger.info(f"Partition {partition_index} has no files to ingest")
            return

        # Create IngestInfo records for each file if backend reporting is enabled
        ingest_file_ids = self._create_ingest_info_records(
            partition_index=partition_index,
            partition_uris=partition_uris,
            omics_dataset_name=omics_dataset_name,
            parent_ingest_id=parent_ingest_id,
        )

        # Report STARTED status for each created IngestInfo record
        self._report_ingest_started(ingest_file_ids=ingest_file_ids)

        # Download files and perform ingest
        first_bucket, _ = self._parse_gcs_uri(partition_uris[0])
        workspace_manager = WorkspaceFileManager(bucket_name=first_bucket)
        with workspace_manager.temp_workspace() as workspace:
            local_paths = self._download_partition_files(
                partition_index=partition_index,
                partition_uris=partition_uris,
                workspace=workspace,
            )

            try:
                self.ingestor.ingest_h5ads_partition(
                    ingest_plan=ingest_plan,
                    local_h5ad_paths=local_paths,
                )
                # Report SUCCEEDED status for all files
                self._report_ingest_succeeded(ingest_file_ids=ingest_file_ids, partition_index=partition_index)
            except Exception as e:
                # Report FAILED status for all files
                self._report_ingest_failed(
                    ingest_file_ids=ingest_file_ids, error_message=str(e), partition_index=partition_index
                )
                raise
