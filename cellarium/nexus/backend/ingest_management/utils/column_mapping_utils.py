from typing import TYPE_CHECKING, Iterable

import pydantic

from cellarium.nexus.omics_datastore.bq_avro_schemas import cell_management

if TYPE_CHECKING:
    from cellarium.nexus.backend.ingest_management import models


def _choices_from_schema(*, model: type[pydantic.BaseModel], exclude: Iterable[str]) -> list[tuple[str, str]]:
    """
    Build choices from a Pydantic sclearchema's fields using their titles.

    :param model: Pydantic model class to introspect
    :param exclude: Iterable of field names to exclude

    :raise ValueError: If field metadata retrieval fails

    :return: List of tuples in the form ``(field_name, title)``
    """
    exclude_set = set(exclude)
    choices: list[tuple[str, str]] = []
    # Pydantic v2 stores fields on `model_fields`
    for name, field in model.model_fields.items():
        if name in exclude_set:
            continue
        title = field.title or name
        choices.append((name, title))
    return choices


def get_obs_column_choices() -> list[tuple[str, str]]:
    """
    Get available column choices for obs mapping from ``CellInfoBQAvroSchema``.

    :raise ValueError: If field metadata retrieval fails

    :return: List of tuples in the form ``(field_name, title)``
    """
    return _choices_from_schema(
        model=cell_management.CellInfoBQAvroSchema,
        exclude=["id", "ingest_id", "metadata_extra", "tag"],
    )


def get_var_column_choices() -> list[tuple[str, str]]:
    """
    Get available column choices for var mapping from ``FeatureInfoBQAvroSchema``.

    :raise ValueError: If field metadata retrieval fails

    :return: List of tuples in the form ``(field_name, title)``
    """
    return _choices_from_schema(
        model=cell_management.FeatureInfoBQAvroSchema,
        exclude=["id", "metadata_extra"],
    )


def create_column_mapping(column_mapping_obj: "models.ColumnMapping | None") -> dict | None:
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
