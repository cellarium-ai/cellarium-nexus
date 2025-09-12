from typing import Iterable

import pydantic
from cellarium.nexus.omics_datastore.bq_avro_schemas import cell_management


def _choices_from_schema(
    *, model: type[pydantic.BaseModel], exclude: Iterable[str]
) -> list[tuple[str, str]]:
    """
    Build choices from a Pydantic schema's fields using their titles.

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
