import tempfile

from cellarium.nexus.clients.nexus_backend_client import NexusBackendAPIClient
from cellarium.nexus.omics_datastore.bq_ops import BigQueryDataValidator
from cellarium.nexus.omics_datastore.bq_ops.validate import validate_remote_file_size
from cellarium.nexus.shared.utils import gcp


class NexusDataValidationCoordinator:
    def __init__(self, *, nexus_backend_api_url: str, validation_report_id: int, max_bytes_valid_per_file: int) -> None:
        """
        Initialize the Nexus data controller.

        :param nexus_backend_api_url: URL for Nexus backend API
        :param validation_report_id: ID of validation report
        """
        self.backend_client = NexusBackendAPIClient(api_url=nexus_backend_api_url)
        self.validation_report_id = validation_report_id
        self.max_bytes_valid_per_file = max_bytes_valid_per_file

    def validate_and_report(self, adata_gcs_path: str, validation_methods: list[str]):
        """
        Validate a single AnnData file and report validation results.

        :param adata_gcs_path: GCS path to the AnnData file to validate
        :param validation_methods: List of validation method names to apply

        :raise: ValidationError if validation fails
        """
        size_is_valid = validate_remote_file_size(
            adata_gcs_path=adata_gcs_path, max_size_bytes=self.max_bytes_valid_per_file
        )
        size_validation_message = "File size is too large." if not size_is_valid else None
        self.backend_client.create_validation_report_item(
            report_id=self.validation_report_id,
            input_file_path=adata_gcs_path,
            validator_name="nexus.omics_datastore.bq_ops.validate.validate_remote_file_size",
            is_valid=size_is_valid,
            message=size_validation_message,
        )
        if not size_is_valid:
            return

        with tempfile.TemporaryDirectory() as temp_dir:
            bucket_name, gcs_file_path = gcp.get_bucket_name_and_file_path_from_gc_path(full_gs_path=adata_gcs_path)
            anndata_file_local_path = f"{temp_dir}/adata.h5ad"

            gcp.download_file_from_bucket(
                bucket_name=bucket_name,
                source_blob_name=gcs_file_path,
                destination_file_name=anndata_file_local_path,
            )

            for validation_method in validation_methods:
                is_valid, error_messages, has_warnings = BigQueryDataValidator.call_validation_method(
                    validation_method=validation_method, anndata_file_local_path=anndata_file_local_path
                )
                message = "\n ".join(error_messages) if error_messages else None
                self.backend_client.create_validation_report_item(
                    report_id=self.validation_report_id,
                    input_file_gcs_path=adata_gcs_path,
                    validator_name=validation_method,
                    is_valid=is_valid,
                    message=message,
                )

    def validate_and_report_multiple(self, adata_gcs_paths: list[str], validation_methods: list[str]):
        """
        Validate multiple AnnData files and report validation results for each.

        :param adata_gcs_paths: List of GCS paths to AnnData files to validate
        :param validation_methods: List of validation method names to apply to each file

        :raise: ValidationError if validation fails for any file
        """
        for adata_gcs_path in adata_gcs_paths:
            self.validate_and_report(adata_gcs_path=adata_gcs_path, validation_methods=validation_methods)
