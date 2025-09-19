from importlib import import_module
from pathlib import Path

from cellarium.nexus.omics_datastore.bq_ops.validate import validate_remote_file_size


def import_from_string(dotted_path: str):
    module_path, func_name = dotted_path.rsplit(".", 1)
    module = import_module(module_path)
    return getattr(module, func_name)


class BigQueryDataValidator:
    @staticmethod
    def validate_remote_file_size(adata_gcs_path: str, max_size_bytes: int, raise_exception=True) -> bool:
        is_valid = validate_remote_file_size(adata_gcs_path=adata_gcs_path, max_size_bytes=max_size_bytes)

        if not is_valid and raise_exception:
            raise ValueError(
                f"File size is too large. Expected <= {max_size_bytes} bytes. You can split input file into multiple "
                f"and try again."
            )

        return is_valid

    @staticmethod
    def call_validation_method(validation_method: str, anndata_file_local_path: str) -> tuple[bool, list[str], bool]:
        validation_method_func = import_from_string(dotted_path=validation_method)
        return validation_method_func(h5ad_path=anndata_file_local_path)

    @staticmethod
    def call_validation_methods(
        validation_methods: list[str], anndata_file_local_path: str | Path
    ) -> list[tuple[bool, list[str], bool]]:
        validation_results = []
        for validation_method in validation_methods:
            validation_method_func = import_from_string(dotted_path=validation_method)
            validation_results.append(validation_method_func(h5ad_path=anndata_file_local_path))

        return validation_results
