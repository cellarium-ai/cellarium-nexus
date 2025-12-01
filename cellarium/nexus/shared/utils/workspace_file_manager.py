"""
Manage temporary workspaces and file operations for data pipelines.
"""

import json
import logging
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from cellarium.nexus.shared.utils import gcp

logger = logging.getLogger(__name__)


class WorkspaceFileManager:
    """
    Manage temporary workspaces and file operations.

    Handle temporary directory creation, cloud storage operations
    """

    def __init__(self, *, bucket_name: str) -> None:
        """
        Initialize the workspace file manager.

        :param bucket_name: Cloud storage bucket name for all operations

        :raise ValueError: If bucket_name is empty
        """
        if not bucket_name:
            raise ValueError("bucket_name cannot be empty")
        self.bucket_name = bucket_name

    @contextmanager
    def temp_workspace(self) -> Iterator[dict[str, Path]]:
        """
        Create a temporary workspace with organized directories.

        Yield a dictionary with 'root', 'input', and 'output' paths.
        Automatically clean up on exit.

        :raise OSError: If directory creation fails

        :return: Context manager yielding dict of paths
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            paths = {
                "root": root,
                "input": root / "input_files",
                "output": root / "output_files",
            }

            paths["input"].mkdir(exist_ok=True)
            paths["output"].mkdir(exist_ok=True)

            logger.info(f"Created temporary workspace at {root}")
            yield paths
            logger.debug(f"Cleaning up temporary workspace at {root}")

    def download_file_from_bucket(
        self,
        *,
        remote_path: str,
        local_path: Path,
    ) -> None:
        """
        Download a file from cloud storage to local path.

        :param remote_path: Path to file within bucket (e.g., "path/to/file.txt")
        :param local_path: Local filesystem path

        :raise IOError: If download fails
        """
        logger.info(f"Downloading gs://{self.bucket_name}/{remote_path} to {local_path}")
        gcp.download_file_from_bucket(
            bucket_name=self.bucket_name,
            source_blob_name=remote_path,
            destination_file_name=local_path,
        )

    def upload_file_to_bucket(
        self,
        *,
        local_path: Path,
        remote_path: str,
    ) -> None:
        """
        Upload a single file to cloud storage.

        :param local_path: Local filesystem path
        :param remote_path: Path within bucket (e.g., "path/to/file.txt")

        :raise IOError: If upload fails
        """
        logger.info(f"Uploading file {local_path} to gs://{self.bucket_name}/{remote_path}")
        gcp.upload_file_to_bucket(
            local_file_path=str(local_path),
            bucket_name=self.bucket_name,
            blob_name=remote_path,
        )

    def upload_directory_to_bucket(
        self,
        *,
        local_path: Path,
        remote_path: str,
        max_workers: int = 8,
    ) -> None:
        """
        Upload an entire directory to cloud storage concurrently in a process pool.

        :param local_path: Local directory path
        :param remote_path: Directory path within bucket (e.g., "path/to/dir")
        :param max_workers: Number of workers for concurrent uploads, defaults to 8.

        :raise IOError: If upload fails
        """
        logger.info(f"Uploading directory {local_path} to gs://{self.bucket_name}/{remote_path}")
        gcp.transfer_directory_to_bucket(
            bucket_name=self.bucket_name,
            local_directory_path=local_path,
            prefix=remote_path,
            max_workers=max_workers,
        )

    def delete_directory_from_bucket(
        self,
        *,
        remote_path: str,
    ) -> list[str]:
        """
        Delete a directory and all files under it from cloud storage.

        :param remote_path: Directory path to delete (e.g., "path/to/dir")

        :raise IOError: If deletion fails

        :return: List of deleted file paths
        """
        logger.info(f"Deleting directory gs://{self.bucket_name}/{remote_path}")
        deleted_files = gcp.delete_files_from_bucket(
            bucket_name=self.bucket_name,
            prefix=remote_path,
        )
        logger.info(f"Deleted {len(deleted_files)} files from gs://{self.bucket_name}/{remote_path}")
        return deleted_files

    def load_json_from_bucket(
        self,
        *,
        remote_path: str,
    ) -> dict[str, Any]:
        """
        Load and deserialize JSON from cloud storage.

        :param remote_path: Path to JSON file within bucket (e.g., "path/to/data.json")

        :raise IOError: If file cannot be read
        :raise ValueError: If file is not valid JSON

        :return: Deserialized JSON data as dictionary
        """
        logger.info(f"Loading JSON from gs://{self.bucket_name}/{remote_path}")
        content = gcp.download_blob_as_string(
            bucket_name=self.bucket_name,
            blob_name=remote_path,
        )
        return json.loads(content)

    def save_json_to_bucket(
        self,
        *,
        data: dict[str, Any],
        remote_path: str,
    ) -> str:
        """
        Serialize and save JSON data to cloud storage.

        :param data: Dictionary to serialize as JSON
        :param remote_path: Destination path within bucket (e.g., "path/to/data.json")

        :raise IOError: If file operations fail

        :return: Full cloud storage path (gs://bucket/path)
        """
        logger.info(f"Saving JSON to gs://{self.bucket_name}/{remote_path}")
        content = json.dumps(data, indent=2)
        gcp.upload_string_as_blob(
            content=content,
            bucket_name=self.bucket_name,
            blob_name=remote_path,
        )
        return f"gs://{self.bucket_name}/{remote_path}"
