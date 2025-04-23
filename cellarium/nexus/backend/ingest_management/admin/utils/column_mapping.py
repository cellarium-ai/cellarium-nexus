from cellarium.nexus.backend.ingest_management import models


def create_column_mapping(column_mapping_obj: models.ColumnMapping | None) -> dict | None:
    """
    Create a column mapping dictionary from a ColumnMapping model instance.

    :param column_mapping_obj: ColumnMapping model instance or None

    :return: Dictionary with obs and var mappings or None if no mappings exist
    """
    if not column_mapping_obj:
        return None

    result = {}

    obs_mapping = {m.input_column: m.schema_column for m in column_mapping_obj.obs_mappings.all()}
    if obs_mapping:
        result["obs_mapping"] = obs_mapping

    var_mapping = {m.input_column: m.schema_column for m in column_mapping_obj.var_mappings.all()}
    if var_mapping:
        result["var_mapping"] = var_mapping

    return result or None
