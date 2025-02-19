import logging
import pathlib
from urllib.parse import urlparse

from google.cloud import storage
from google.cloud.storage import transfer_manager

logger = logging.getLogger(__name__)


def get_bucket_name_and_file_path_from_gc_path(full_gs_path: str) -> tuple[str, str]:
    parsed = urlparse(full_gs_path)
    bucket_name = parsed.netloc

    # Construct full path including the fragment if it exists
    file_path = parsed.path.lstrip("/")
    if parsed.fragment:
        file_path = f"{file_path}#{parsed.fragment}"

    return bucket_name, file_path


def download_file_from_bucket(bucket_name: str, source_blob_name: str, destination_file_name: pathlib.Path) -> None:
    """
    Download a file from GCS

    :param bucket_name: Bucket name in Google Cloud Storage
    :param source_blob_name: Name of the source blob in Google Cloud Storage Bucket
    :param destination_file_name: Local file name where to save the blob
    """
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)

    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)


def upload_file_to_bucket(local_file_path: str, bucket_name: str, blob_name: str) -> None:
    client = storage.Client()
    bucket_name = client.get_bucket(bucket_or_name=bucket_name)
    blob = bucket_name.blob(blob_name=blob_name)
    blob.upload_from_filename(filename=local_file_path)


def list_blobs(bucket_name: str, prefix: str | None = None) -> list[storage.Blob]:
    """
    List all the blobs in the bucket

    :param bucket_name: Bucket name in Google Cloud Storage
    :param prefix: list all blobs with a prefix (a.k.a. "directory")

    :return: List of blobs with the provided prefix
    """
    storage_client = storage.Client()
    return storage_client.list_blobs(bucket_name, prefix=prefix)


def upload_many_blobs_with_transfer_manager(bucket_name, file_paths, prefix, workers=8):
    """Upload every file in a list to a bucket, concurrently in a process pool.

    Each blob name is derived from the filename, not including the
    `source_directory` parameter. For complete control of the blob name for each
    file (and other aspects of individual blob metadata), use
    transfer_manager.upload_many() instead.
    """
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)

    filenames_blob_pairs = []

    for file_path in file_paths:
        file_name = file_path.split("/")[-1]
        blob = storage.Blob(bucket=bucket, name=f"{prefix}/{file_name}")
        filenames_blob_pairs.append((file_path, blob))

    results = transfer_manager.upload_many(
        file_blob_pairs=filenames_blob_pairs, max_workers=workers, upload_kwargs={"timeout": 3000}
    )

    number_of_uploaded = 0
    for name, result in zip(file_paths, results):
        # The results list is either `None` or an exception for each filename in
        # the input list, in order.

        if isinstance(result, Exception):
            logger.info(f"Failed to upload {name} due to exception: {result}")
        else:
            number_of_uploaded += 1

    logger.info(f"Successfully uploaded {number_of_uploaded} files")


def transfer_directory_to_bucket(
    bucket_name: str, local_directory_path: pathlib.Path, prefix: str, max_workers: int = 8
):
    directory = pathlib.Path(local_directory_path)
    file_paths = [str(file) for file in directory.rglob("*") if file.is_file()]
    upload_many_blobs_with_transfer_manager(
        bucket_name=bucket_name, file_paths=file_paths, prefix=prefix, workers=max_workers
    )
