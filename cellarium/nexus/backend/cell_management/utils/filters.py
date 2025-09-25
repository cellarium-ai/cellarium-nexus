"""
Utility helpers for admin filters.

Provide normalization and validation helpers for filter statements used by the
Cell Info admin endpoints.
"""

from __future__ import annotations

import decimal as py_decimal
import types as py_types
import typing as t


def normalize_filter_statements(*, filter_statements: dict[str, t.Any]) -> dict[str, t.Any]:
    """
    Normalize filter statements for BigQuery where-clause rendering.

    Convert equality operators with list values to membership operators, so the
    SQL helper receives types it supports.

    :param filter_statements: Mapping of "column__op" to values.

    :raise: None

    :return: Normalized mapping of "column__op" to values.
    """
    normalized: dict[str, t.Any] = {}
    for key, value in filter_statements.items():
        if "__" in key:
            col, op = key.rsplit("__", 1)
        else:
            col, op = key, "eq"
        if op in ("eq", "not_eq") and isinstance(value, list):
            if len(value) == 1:
                normalized[f"{col}__{op}"] = value[0]
            else:
                mapped = "in" if op == "eq" else "not_in"
                normalized[f"{col}__{mapped}"] = value
        else:
            normalized[key] = value
    return normalized


def resolve_base_type(anno: t.Any) -> t.Any:
    """
    Resolve typing annotation to a base Python type.

    Handle ``Annotated[T, ...]`` and optional/union types (``Union[T, None]`` or
    ``T | None``). Prefer ``bool`` over ``str`` over numeric when multiple types
    remain in a union.

    :param anno: The typing annotation to resolve.

    :raise: None

    :return: The resolved base type (for example, :class:`str`, :class:`int`, :class:`bool`)
        or the original annotation if unresolved.
    """
    current = anno
    while True:
        origin = t.get_origin(current)
        args = t.get_args(current)
        if origin is None:
            return current
        # Unwrap Annotated[T, ...]
        if origin is t.Annotated and args:
            current = args[0]
            continue
        # Unwrap Optional/Union removing NoneType; if multiple remain, choose a representative
        if origin in (t.Union, py_types.UnionType) or str(origin) in {"typing.Union"}:
            non_none = tuple(a for a in args if a is not type(None))  # noqa: E721
            if len(non_none) == 1:
                current = non_none[0]
                continue
            # Heuristics: prefer bool, then str, then numeric types
            if any(x is bool for x in non_none):
                return bool
            if any(x is str for x in non_none):
                return str
            if any(x in (int, float, py_decimal.Decimal) for x in non_none):
                # Prefer float for numeric unions
                return float
            # Fallback: return as-is (unsupported complex union)
            return current
        return current
