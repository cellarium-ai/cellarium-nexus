"""
Control and manage SOMA extract operations for Nexus.
"""

import logging

from cellarium.nexus.clients import NexusBackendAPIClient
from cellarium.nexus.omics_datastore.soma_ops import TileDBSOMADataOperator
from cellarium.nexus.shared.utils.workspace_file_manager import WorkspaceFileManager

logger = logging.getLogger(__name__)


class SomaDataOpsCoordinator:
    """
    Control and manage SOMA extract operations for Nexus.

    Provide a high-level interface for SOMA extract operations,
    delegating low-level operations to TileDBSOMADataOperator.
    """

    def __init__(self, *, experiment_uri: str, nexus_backend_api_url: str) -> None:
        """
        Initialize the SOMA data operations coordinator.

        :param experiment_uri: URI of the SOMA experiment
        :param nexus_backend_api_url: URL for Nexus backend API

        :raise ValueError: If experiment_uri is empty
        """
        if not experiment_uri:
            raise ValueError("experiment_uri cannot be empty")

        self.experiment_uri = experiment_uri
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url)
        self.soma_operator = TileDBSOMADataOperator(experiment_uri=experiment_uri)

        logger.info(f"Initialized SomaDataOpsCoordinator for {experiment_uri}")

    def _workspace(self, *, bucket_name: str) -> WorkspaceFileManager:
        """
        Create a workspace file manager for the given bucket.

        :param bucket_name: GCS bucket name

        :return: WorkspaceFileManager instance
        """
        return WorkspaceFileManager(bucket_name=bucket_name)
