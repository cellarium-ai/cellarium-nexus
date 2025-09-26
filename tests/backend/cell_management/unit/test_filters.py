import decimal
import typing

from cellarium.nexus.backend.cell_management.utils import filters


def test_normalize_filter_eq_list_to_in_and_singleton_collapses() -> None:
    """
    Normalize equality lists to membership and collapse singletons to scalars.
    """
    data = {
        "organism__eq": ["human", "mouse"],
        "tissue__eq": ["lung"],
        "age__eq": 42,
    }
    out = filters.normalize_filter_statements(filter_statements=data)
    assert out["organism__in"] == ["human", "mouse"]
    assert "organism__eq" not in out
    assert out["tissue__eq"] == "lung"
    assert out["age__eq"] == 42


def test_normalize_filter_not_eq_list_to_not_in() -> None:
    """
    Normalize not-equal lists to not-in while preserving singletons as scalars.
    """
    data = {
        "organism__not_eq": ["human", "mouse"],
        "tissue__not_eq": ["lung"],
    }
    out = filters.normalize_filter_statements(filter_statements=data)
    assert out["organism__not_in"] == ["human", "mouse"]
    assert "organism__not_eq" not in out
    assert out["tissue__not_eq"] == "lung"


def test_normalize_filter_passthrough_for_other_ops() -> None:
    """
    Pass through non-equality operators and plain columns unchanged.
    """
    data = {
        "age__gt": 10,
        "height__lt": 15.5,
        "weight": 70,
    }
    out = filters.normalize_filter_statements(filter_statements=data)
    assert out == data


# resolve_base_type tests


def test_resolve_base_type_unwraps_annotated_and_optional() -> None:
    """
    Unwrap Annotated and Optional to their underlying base primitive type.
    """
    T = typing.Annotated[int, "meta"]
    anno = typing.Optional[T]
    assert filters.resolve_base_type(anno) is int


def test_resolve_base_type_union_precedence_bool_over_str_over_numeric() -> None:
    """
    Prefer bool over str over numeric for unions with multiple primitive types.
    """
    assert filters.resolve_base_type(typing.Union[bool, str]) is bool
    assert filters.resolve_base_type(typing.Union[str, int]) is str
    assert filters.resolve_base_type(typing.Union[int, float, decimal.Decimal]) is float


def test_resolve_base_type_simple_primitives_and_passthrough() -> None:
    """
    Return primitives as-is and pass through complex unions unchanged.
    """
    assert filters.resolve_base_type(int) is int
    assert filters.resolve_base_type(str) is str
    # Unresolvable complex union returns original
    U = typing.Union[tuple, set]
    assert filters.resolve_base_type(U) == U
