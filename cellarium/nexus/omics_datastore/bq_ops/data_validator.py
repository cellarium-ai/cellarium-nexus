from importlib import import_module


def import_from_string(dotted_path: str):
    module_path, func_name = dotted_path.rsplit(".", 1)
    module = import_module(module_path)
    return getattr(module, func_name)


class BigQueryDataValidator:
    @staticmethod
    def call_validation_method(validation_method: str, anndata_file_local_path: str) -> tuple[bool, list[str], bool]:
        validation_method_func = import_from_string(dotted_path=validation_method)
        return validation_method_func(h5ad_path=anndata_file_local_path)
