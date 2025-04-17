from typing import List, Tuple

from django.apps import apps


def get_obs_column_choices() -> List[Tuple[str, str]]:
    """
    Get available column choices for obs mapping from CellInfo model.

    :return: List of (field_name, verbose_name) tuples

    :raise LookupError: If model is not found
    :raise AttributeError: If field metadata is missing
    """
    CellInfo = apps.get_model("cell_management", "CellInfo")
    return [
        (field.name, field.verbose_name)
        for field in CellInfo._meta.fields
        if field.name not in ["id", "ingest", "metadata_extra", "tag"]
    ]


def get_var_column_choices() -> List[Tuple[str, str]]:
    """
    Get available column choices for var mapping from FeatureInfo model.

    :return: List of (field_name, verbose_name) tuples

    :raise LookupError: If model is not found
    :raise AttributeError: If field metadata is missing
    """
    FeatureInfo = apps.get_model("cell_management", "CellFeatureInfo")
    return [
        (field.name, field.verbose_name)
        for field in FeatureInfo._meta.fields
        if field.name not in ["id", "metadata_extra"]
    ]
