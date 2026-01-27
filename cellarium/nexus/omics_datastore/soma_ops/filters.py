"""
Filter translation utilities for SOMA operations.

This module provides functions to convert Nexus filter dicts to SOMA value_filter expressions.
"""

import re

from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaFilterError


def _escape_string_literal(value: str) -> str:
    """
    Escape a string value for use in SOMA value_filter expressions.

    :param value: String value to escape

    :return: Escaped and quoted string literal
    """
    # Escape backslashes and quotes
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _format_literal(value: object) -> str:
    """
    Format a Python value as a SOMA expression literal.

    :param value: Python value to format

    :raise SomaFilterError: If value type is not supported

    :return: SOMA literal representation
    """
    if isinstance(value, str):
        return _escape_string_literal(value)
    elif isinstance(value, bool):
        # Boolean must be checked before int since bool is a subclass of int
        return str(value)
    elif isinstance(value, (int, float)):
        return str(value)
    elif isinstance(value, list):
        formatted_items = [_format_literal(item) for item in value]
        return f"[{', '.join(formatted_items)}]"
    else:
        raise SomaFilterError(f"Unsupported value type for SOMA literal: {type(value).__name__}")


def _parse_filter_key(key: str) -> tuple[str, str]:
    """
    Parse a filter key into column name and operator.

    :param key: Filter key in format 'column__op' or 'alias.column__op'

    :raise SomaFilterError: If key format is invalid

    :return: Tuple of (column_name, operator)
    """
    # Remove table alias prefix if present (e.g., 'c.', 't.')
    key_without_alias = re.sub(r"^[a-zA-Z_]\.", "", key)

    # Split on the last occurrence of '__'
    if "__" not in key_without_alias:
        raise SomaFilterError(f"Invalid filter key format: {key}. Expected 'column__op' format.")

    parts = key_without_alias.rsplit("__", 1)
    if len(parts) != 2:
        raise SomaFilterError(f"Invalid filter key format: {key}. Expected 'column__op' format.")

    column, operator = parts
    if not column or not operator:
        raise SomaFilterError(f"Invalid filter key format: {key}. Column and operator cannot be empty.")

    return column, operator


def _build_condition(column: str, operator: str, value: object) -> str:
    """
    Build a single SOMA filter condition.

    :param column: Column name
    :param operator: Operator suffix (e.g., 'eq', 'in', 'gt')
    :param value: Filter value

    :raise SomaFilterError: If operator is not supported

    :return: SOMA filter condition expression
    """
    if operator == "eq":
        return f"({column} == {_format_literal(value)})"
    elif operator == "not_eq":
        return f"({column} != {_format_literal(value)})"
    elif operator == "in":
        if not isinstance(value, list):
            raise SomaFilterError(f"Operator 'in' requires a list value, got {type(value).__name__}")
        return f"({column} in {_format_literal(value)})"
    elif operator == "not_in":
        if not isinstance(value, list):
            raise SomaFilterError(f"Operator 'not_in' requires a list value, got {type(value).__name__}")
        return f"({column} not in {_format_literal(value)})"
    elif operator == "gt":
        return f"({column} > {_format_literal(value)})"
    elif operator == "gte":
        return f"({column} >= {_format_literal(value)})"
    elif operator == "lt":
        return f"({column} < {_format_literal(value)})"
    elif operator == "lte":
        return f"({column} <= {_format_literal(value)})"
    else:
        raise SomaFilterError(f"Unsupported operator: {operator}")


def build_soma_value_filter(filters: dict[str, object] | None) -> str:
    """
    Build a SOMA value_filter expression from a filter dict.

    Convert Nexus filter dict format to SOMA value_filter expression string.
    Supported operators: eq, not_eq, in, not_in, gt, gte, lt, lte.

    :param filters: Dict of filter conditions using the Nexus format

    :raise SomaFilterError: If an unsupported operator is used or filter format is invalid

    :return: SOMA value_filter expression string
    """
    if not filters:
        return ""

    conditions = []
    for key, value in filters.items():
        column, operator = _parse_filter_key(key)
        condition = _build_condition(column, operator, value)
        conditions.append(condition)

    # Combine all conditions with 'and'
    if len(conditions) == 1:
        return conditions[0]
    else:
        return " and ".join(conditions)
