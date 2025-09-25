"""
Utility functions for checking curriculum existence.
"""

from django.conf import settings
from google.cloud import storage


def check_curriculum_exists(name: str) -> bool:
    """
    Check if a curriculum with the given name already exists in the bucket.

    This function checks if there are any files in the data-extracts directory
    with the provided curriculum name, which would indicate that a curriculum
    with this name already exists.

    :param name: Name of the curriculum to check

    :raise: google.cloud.exceptions.GoogleCloudError: If there's an error communicating with Google Cloud Storage

    :return: True if curriculum exists, False otherwise
    """
    prefix = f"{settings.BACKEND_PIPELINE_DIR}/{settings.PIPELINE_DATA_EXTRACTS_DIR}/{name}"

    storage_client = storage.Client()
    blobs = list(storage_client.bucket(settings.BUCKET_NAME_PRIVATE).list_blobs(prefix=prefix, max_results=1))

    return len(blobs) > 0
