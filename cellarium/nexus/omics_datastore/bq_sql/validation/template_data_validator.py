import typing as t

from nexus.omics_datastore.bq_sql.constants import ComparisonOperators
from nexus.omics_datastore.bq_sql.validation import exceptions


def validate_not_empty(value: str | list | dict) -> None:
    """
     Validate that a given `value` is not empty.

    :param value: The value to be validated.
    :raises exceptions.TemplateDataValidationError: Raised if `value` is empty.
    """
    if not value:
        raise exceptions.TemplateDataValidationError("Value should not be an empty")


def validate_sql_filter(filter_name: str, filter_value: list | str | bool | int | float) -> None:
    """
    Validates the SQL filter parameters based on the type of filter applied, ensuring adherence to the supported
    filter types and value constraints.

    This function checks the 'filter_name' to ascertain if the filter is of type 'eq' (exact) or 'in' ('in' comparison).
    It then validates 'filter_value' based on the determined filter type. For 'eq', 'filter_value' must be a basic data
    type (int, str, float, or bool), while for 'in', 'filter_value' must be a list containing elements of the same type.

    :param filter_name: The name of the filter, expected in the '<column_name>__<filter_type>' format
    :param filter_value: The value used for filtering
    :raises exceptions.TemplateDataValidationError: If validation fails
    """
    column_name, filter_type = filter_name.split("__")
    if filter_type not in ComparisonOperators.CURRENTLY_SUPPORTED:
        raise exceptions.TemplateDataValidationError(
            f"`Only <column_name>__eq or <column_name>__in available as types of filters, got {filter_type}"
        )

    if filter_type == ComparisonOperators.EQUAL:
        if not isinstance(filter_value, (int, str, float, bool)):
            raise exceptions.TemplateDataValidationError(
                f"When filtering by the exact value, `filter_value` supposed to be a python basic "
                f"data type in {filter_name} filter"
            )
    elif filter_type == ComparisonOperators.IN:
        if not isinstance(filter_value, list):
            raise exceptions.TemplateDataValidationError(
                f"When filtering by the range of values, `filter_value` supposed to be a list with elements of the same "
                f"type in {filter_name} filter"
            )

        types = set(map(type, filter_value))
        is_single_type_list = len(types) == 1

        if not is_single_type_list:
            raise exceptions.TemplateDataValidationError(
                f"All elements in the list must be of the same type for {filter_name} filter"
            )

        filter_values_type = types.pop()

        if filter_values_type not in [str, int, float]:
            raise exceptions.TemplateDataValidationError(
                f"Only int, str, and floats allowed for interval filter, got {filter_values_type}"
            )


def validate_filters(filters: dict[str, str]) -> None:
    """
    Validate the structure and contents of filter statements provided for SQL query generation.

    :param filters: A dictionary representing filter statements for the SQL query.
    :raises exceptions.TemplateDataValidationError: If any filter statement doesn't pass the validation rules.
    """
    for filter_name, filter_value in filters.items():
        validate_sql_filter(filter_name=filter_name, filter_value=filter_value)


def validate_column_name(column_name: str) -> None:
    """
    Validate column name. Check for number of periods in column name as well as for an empty string

    :param column_name: The column name to validate.
    :raises ValueError: If the column name contains more than one period or is an empty string.
    """
    column_split = column_name.split(".")

    if len(column_split) > 2:
        raise ValueError(
            f"Column {column_name} has more than one period in its name. It can contain only one period, "
            f"which separates the alias table name from the column name itself."
        )

    if not column_name:
        raise ValueError("Empty strings are not allowed for column names")
