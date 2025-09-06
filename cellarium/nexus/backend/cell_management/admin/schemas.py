"""
Define Pydantic schemas for Cell Info admin endpoints.
"""

from __future__ import annotations

import typing as t

import pydantic as pydantic


class FilterItem(pydantic.BaseModel):
    """
    Describe a single filter condition.

    :param field: Field key the filter applies to
    :param operator: Operator to apply (eq, in, gt, etc.)
    :param value: Filter value; for categorical provide list[str], for numbers provide float, for booleans provide bool

    :raise: pydantic.ValidationError

    :return: None
    """

    field: str
    operator: str
    value: t.Union[list[str], float, bool]


class FiltersPayload(pydantic.BaseModel):
    """
    Represent the filters request payload sent from the front-end.

    :param dataset: Selected dataset name
    :param filters: Filter conditions; accept either a mapping of
        ``{"column__op": value}`` or a list of FilterItem objects.

    :raise: pydantic.ValidationError

    :return: None
    """

    dataset: str
    filters: t.Union[dict[str, t.Any], list[FilterItem]] = pydantic.Field(default_factory=dict)

    def to_filter_statements(self) -> dict[str, t.Any]:
        """
        Convert filters to a BigQuery-compatible statements mapping.

        :raise: pydantic.ValidationError

        :return: A dictionary with keys in the form ``column__op`` mapping to values.
        """
        if isinstance(self.filters, dict):
            return self.filters
        statements: dict[str, t.Any] = {}
        for item in self.filters:
            key = f"{item.field}__{item.operator}"
            statements[key] = item.value
        return statements


class CountResponse(pydantic.BaseModel):
    """
    Represent the count endpoint response.

    :param count: Mocked count value

    :raise: pydantic.ValidationError

    :return: None
    """

    count: int


class FieldMeta(pydantic.BaseModel):
    """
    Describe a filterable field and UI hints.

    :param key: Internal field key
    :param label: Human-readable label
    :param type: Field type; allowed: string, number, boolean, text
    :param operators: Allowed operators
    :param suggest: Whether the field supports suggestions
    :param suggest_mode: Suggestion matching mode (for example, prefix)
    :param suggest_min_chars: Minimum characters to trigger suggestions

    :raise: pydantic.ValidationError

    :return: None
    """

    key: str
    label: str
    type: str
    operators: list[str] = pydantic.Field(default_factory=list)
    suggest: bool = False
    suggest_mode: t.Optional[str] = None
    suggest_min_chars: t.Optional[int] = None
