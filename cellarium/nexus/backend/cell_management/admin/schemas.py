"""
Define Pydantic schemas for Cell Info admin endpoints.
"""

from __future__ import annotations

import typing as t

import pydantic as pydantic


class FiltersPayload(pydantic.BaseModel):
    """
    Represent the filters request payload sent from the front-end.

    :param dataset: Selected dataset name
    :param filters: Filter conditions as a mapping of ``{"column__op": value}``.

    :raise: pydantic.ValidationError

    :return: None
    """

    dataset: str
    filters: dict[str, t.Any] = pydantic.Field(default_factory=dict)


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
