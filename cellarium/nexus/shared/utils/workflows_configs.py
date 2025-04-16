"""
Utility functions for handling workflow configurations.
"""

import datetime
import os
import tempfile
from typing import TypeVar
import secrets

import yaml
from pydantic import BaseModel
from smart_open import open

from cellarium.nexus.shared import utils

T = TypeVar("T", bound=BaseModel)


def read_component_config(gcs_path: str, schema_class: type[T]) -> T:
    """
    Read a YAML configuration file for a Kubeflow component and convert it to a Pydantic model.

    :param gcs_path: Full GCS path to the YAML file (e.g., 'gs://bucket-name/path/to/file.yaml').
    :param schema_class: Pydantic model class that defines the expected schema.

    :raise: ValueError - If the GCS path is invalid or empty.
    :raise: FileNotFoundError - If the file doesn't exist in the bucket.
    :raise: yaml.YAMLError - If the YAML file is malformed.
    :raise: pydantic.ValidationError - If the YAML data doesn't match the expected schema.

    :return: An instance of the provided schema_class containing the validated configuration.
    """
    with open(gcs_path, "r") as f:
        data = yaml.safe_load(f)
        return schema_class(**data)


def dump_configs_to_bucket(configs: list[BaseModel], bucket_path: str, workers: int = 8) -> list[str]:
    """
    Dump a list of Pydantic model instances to a temporary directory and upload them to a GCS bucket.

    Each file is named using the model's class name (lowercase) and a timestamp.
    Files are uploaded concurrently using transfer_manager for better performance.

    :param configs: Pydantic model instances to dump.
    :param bucket_path: GCS path where the configs should be uploaded (e.g., 'gs://bucket-name/path/to/configs/').
    :param workers: Number of workers for concurrent uploads, defaults to 8.

    :raise: ValueError - If the bucket path is invalid or empty.
    :raise: IOError - If there's an error writing to the temporary directory.
    :raise: Exception - If there's an error uploading to GCS.

    :return: List of GCS paths where the configs were uploaded.
    """
    if not configs:
        return []

    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

    bucket_name, prefix = utils.gcp.get_bucket_name_and_file_path_from_gc_path(bucket_path)

    if prefix:
        # Remove any leading/trailing slashes and add a single trailing slash
        prefix = prefix.strip("/") + "/"

    with tempfile.TemporaryDirectory() as temp_dir:
        local_file_paths = []
        blob_names = []

        for config in configs:
            model_name = config.__class__.__name__.lower()

            filename = f"{model_name}_{timestamp}_{secrets.token_hex(6)}.yaml"
            local_path = os.path.join(temp_dir, filename)

            with open(local_path, "w") as f:
                yaml.dump(config.model_dump(), f)

            local_file_paths.append(local_path)
            blob_names.append(f"{prefix}{filename}")

        utils.gcp.upload_many_blobs_with_transfer_manager(
            bucket_name=bucket_name, file_paths=local_file_paths, prefix=prefix, workers=workers
        )

        uploaded_paths = [f"gs://{bucket_name}/{blob_name}" for blob_name in blob_names]

        return uploaded_paths
