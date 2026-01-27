import pytest

from cellarium.nexus.omics_datastore.soma_ops import SomaFilterError
from cellarium.nexus.omics_datastore.soma_ops import filters as filters_module


def test_build_soma_value_filter_empty_filters() -> None:
    """
    Verify that None and empty dict return empty string.
    """
    assert filters_module.build_soma_value_filter(filters=None) == ""
    assert filters_module.build_soma_value_filter(filters={}) == ""


@pytest.mark.parametrize(
    "filters,expected",
    [
        ({"organism__eq": "Homo sapiens"}, '(organism == "Homo sapiens")'),
        ({"organism__not_eq": "Mus musculus"}, '(organism != "Mus musculus")'),
    ],
)
def test_build_soma_value_filter_eq_ne_operators(filters: dict[str, object], expected: str) -> None:
    """
    Verify __eq and __not_eq operators translate correctly.
    """
    result = filters_module.build_soma_value_filter(filters=filters)
    assert result == expected


@pytest.mark.parametrize(
    "filters,expected",
    [
        ({"age__gt": 30}, "(age > 30)"),
        ({"age__gte": 30}, "(age >= 30)"),
        ({"age__lt": 50}, "(age < 50)"),
        ({"age__lte": 50}, "(age <= 50)"),
    ],
)
def test_build_soma_value_filter_numeric_operators(filters: dict[str, object], expected: str) -> None:
    """
    Verify numeric comparison operators (gt, gte, lt, lte).
    """
    result = filters_module.build_soma_value_filter(filters=filters)
    assert result == expected


@pytest.mark.parametrize(
    "filters,expected",
    [
        ({"is_primary__eq": True}, "(is_primary == True)"),
        ({"is_primary__eq": False}, "(is_primary == False)"),
    ],
)
def test_build_soma_value_filter_boolean_values(filters: dict[str, object], expected: str) -> None:
    """
    Verify boolean values translate to True/False literals.
    """
    result = filters_module.build_soma_value_filter(filters=filters)
    assert result == expected


@pytest.mark.parametrize(
    "filters,expected",
    [
        ({"cell_type__in": ["T Cell", "B Cell"]}, '(cell_type in ["T Cell", "B Cell"])'),
        ({"cell_type__not_in": ["T Cell", "B Cell"]}, '(cell_type not in ["T Cell", "B Cell"])'),
        ({"donor_id__in": [1, 2, 3]}, "(donor_id in [1, 2, 3])"),
    ],
)
def test_build_soma_value_filter_in_not_in_operators(filters: dict[str, object], expected: str) -> None:
    """
    Verify __in and __not_in operators with list values.
    """
    result = filters_module.build_soma_value_filter(filters=filters)
    assert result == expected


def test_build_soma_value_filter_multiple_conditions() -> None:
    """
    Verify multiple filter conditions combine with 'and'.
    """
    result = filters_module.build_soma_value_filter(
        filters={
            "tissue__eq": "lung",
            "cell_type__in": ["T Cell", "B Cell"],
        }
    )
    # Should combine with 'and'
    assert '(tissue == "lung")' in result
    assert '(cell_type in ["T Cell", "B Cell"])' in result
    assert " and " in result


@pytest.mark.parametrize(
    "filters,expected",
    [
        ({"c.organism__eq": "Homo sapiens"}, '(organism == "Homo sapiens")'),
        ({"t.cell_type__in": ["T", "B"]}, '(cell_type in ["T", "B"])'),
    ],
)
def test_build_soma_value_filter_strips_table_alias(filters: dict[str, object], expected: str) -> None:
    """
    Verify table alias prefixes (e.g., 'c.', 't.') are stripped from column names.
    """
    result = filters_module.build_soma_value_filter(filters=filters)
    assert result == expected


def test_build_soma_value_filter_string_escaping() -> None:
    """
    Verify string values with quotes and backslashes are properly escaped.
    """
    result = filters_module.build_soma_value_filter(filters={"tissue__eq": 'lung "left"\\right'})
    # Should escape quotes and backslashes
    assert "lung" in result
    assert "left" in result
    assert "right" in result


def test_build_soma_value_filter_unsupported_operator() -> None:
    """
    Verify unsupported operator raises SomaFilterError.
    """
    with pytest.raises(SomaFilterError):
        filters_module.build_soma_value_filter(filters={"age__between": [1, 2]})


def test_build_soma_value_filter_missing_operator() -> None:
    """
    Verify missing operator suffix raises SomaFilterError.
    """
    with pytest.raises(SomaFilterError):
        filters_module.build_soma_value_filter(filters={"organism": "human"})


@pytest.mark.parametrize(
    "filters",
    [
        {"__eq": "value"},
        {"organism__": "value"},
    ],
)
def test_build_soma_value_filter_empty_column_or_operator(filters: dict[str, object]) -> None:
    """
    Verify empty column name or operator raises SomaFilterError.
    """
    with pytest.raises(SomaFilterError):
        filters_module.build_soma_value_filter(filters=filters)


@pytest.mark.parametrize(
    "filters",
    [
        {"cell_type__in": "T Cell"},
        {"cell_type__not_in": "T Cell"},
    ],
)
def test_build_soma_value_filter_in_not_in_with_non_list(filters: dict[str, object]) -> None:
    """
    Verify __in and __not_in operators with non-list value raise SomaFilterError.
    """
    with pytest.raises(SomaFilterError):
        filters_module.build_soma_value_filter(filters=filters)


def test_build_soma_value_filter_unsupported_value_type() -> None:
    """
    Verify unsupported value type (e.g., set) raises SomaFilterError.
    """
    with pytest.raises(SomaFilterError):
        filters_module.build_soma_value_filter(filters={"cell_type__eq": {"T", "B"}})
