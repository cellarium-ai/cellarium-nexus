import os

import anndata
from cellarium_schema.validate import validate
from smart_open import open

from cellarium.nexus.shared.utils import gcp


def cellarium_validate_gencode_43(h5ad_path: str):
    return validate(h5ad_path, ignore_labels=True, gencode_version=43)


def cellarium_validate_gencode_44(h5ad_path: str):
    return validate(h5ad_path, ignore_labels=True, gencode_version=44)


def validate_remote_file_size(adata_gcs_path: str, max_size_bytes: int) -> bool:
    """
    Validate that an AnnData file does not exceed the specified size limit.

    :param adata_gcs_path: GCS path to the AnnData file (e.g., 'gs://bucket-name/path/to/file.h5ad').
    :param max_size_bytes: Maximum allowed file size in bytes.

    :raise: ValueError - If the GCS path is invalid or empty.
    :raise: google.cloud.exceptions.GoogleCloudError - If there's an error communicating with Google Cloud Storage.
    :raise: google.cloud.exceptions.NotFound - If the specified file doesn't exist.

    :return: True if the file size is within the limit, False otherwise.
    """
    bucket_name, file_path = gcp.get_bucket_name_and_file_path_from_gc_path(full_gs_path=adata_gcs_path)
    file_size = gcp.get_blob_size(bucket_name=bucket_name, file_path=file_path)

    return file_size <= max_size_bytes


def validate_raw_counts(
    h5ad_path: str | bytes | os.PathLike,
    add_labels_file: str = None,
    gencode_version: int = 44,
    ignore_labels: bool = False,
) -> tuple[bool, list, bool]:
    """
    Validate that the AnnData object's .X matrix contains raw counts (not normalized values).

    The function checks if the sum of the matrix is a whole number or very close to a whole number,
    which is typically the case for raw count data but not for normalized data.

    :param h5ad_path: Path to the h5ad file to validate.
    :param add_labels_file: Optional path to additional labels file, not used in this validation.
    :param gencode_version: GENCODE version, not used in this validation.
    :param ignore_labels: Whether to ignore labels, not used in this validation.

    :raise: ValueError - If the AnnData object's .X matrix is empty or None.
    :raise: FileNotFoundError - If the h5ad file cannot be found.

    :return: Tuple containing (is_valid, error_messages, has_warnings).
             is_valid: True if the matrix appears to contain raw counts.
             error_messages: List of error messages if validation fails.
             has_warnings: True if there are warnings but validation passes.
    """
    errors = []
    warnings = False
    tolerance = 0.01

    try:
        # Read the h5ad file
        with open(h5ad_path, "rb") as f:
            adata = anndata.read_h5ad(f)

        # Check if X matrix exists and is not empty
        if adata.X is None or adata.X.size == 0:
            errors.append("AnnData object's .X matrix is empty or None")
            return False, errors, warnings

        # Calculate the sum of all values in the matrix
        matrix_sum = adata.X.sum()

        # Check if the sum is close to an integer
        is_close_to_integer = abs(matrix_sum - round(matrix_sum)) <= tolerance

        if not is_close_to_integer:
            errors.append("The matrix does not appear to contain raw counts. The sum is not close to an integer.")
            return False, errors, warnings

        return True, [], warnings

    except FileNotFoundError:
        errors.append(f"File not found: {h5ad_path}")
        return False, errors, warnings
    except Exception as e:
        errors.append(f"Error validating raw counts: {str(e)}")
        return False, errors, warnings
