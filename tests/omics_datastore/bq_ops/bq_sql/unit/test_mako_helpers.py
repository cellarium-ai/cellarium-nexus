import pytest

from cellarium.nexus.omics_datastore.bq_ops.bq_sql import exceptions as bq_exceptions
from cellarium.nexus.omics_datastore.bq_ops.bq_sql import mako_helpers


def test_build_concat_expression_basic() -> None:
    """
    Build a concat expression with default delimiter and alias.
    """
    expr = mako_helpers.build_concat_expression(columns=["c.a", "c.b"], alias="ab")
    assert expr == "concat(c.a, '##', c.b) as ab"


def test_build_concat_expression_custom_delim() -> None:
    """
    Build a concat expression using a custom delimiter across three columns.
    """
    expr = mako_helpers.build_concat_expression(columns=["c.a", "c.b", "c.c"], alias="abc", delim="|")
    assert expr == "concat(c.a, '|', c.b, '|', c.c) as abc"


def test_build_concat_expression_empty_columns_raises() -> None:
    """
    Raise an error when attempting to build a concat expression with no columns.
    """
    with pytest.raises(ValueError):
        mako_helpers.build_concat_expression(columns=[], alias="x")


def test_where_supported_filters_formatting() -> None:
    """
    Render a where clause covering all supported operators with proper formatting.
    """
    filters = {
        "organism__eq": "Homo sapiens",
        "id__in": [1, 3],
        "flag__not_eq": True,
        "count__gt": 10,
        "ratio__gte": 1.5,
        "score__lt": 0,
        "size__lte": 5,
    }

    result = mako_helpers.where(filters=filters)

    assert result == (
        "where\n"
        "    organism = 'Homo sapiens'\n"
        "    and id in (1, 3)\n"
        "    and flag != TRUE\n"
        "    and count > 10\n"
        "    and ratio >= 1.5\n"
        "    and score < 0\n"
        "    and size <= 5"
    )


def test_where_in_with_strings_quotes() -> None:
    """
    Quote string values for an IN list within the where clause.
    """
    result = mako_helpers.where(filters={"cell_type__in": ["T cell", "neuron"]})
    assert result == "where\n    cell_type in ('T cell', 'neuron')"


def test_where_in_with_floats_now_supported() -> None:
    """
    Render a where clause with an IN list of floats.
    """
    result = mako_helpers.where(filters={"ratio__in": [0.1, 1.5]})
    assert result == "where\n    ratio in (0.1, 1.5)"


def test_where_unsupported_filter_type_raises() -> None:
    """
    Raise an error for an unsupported filter operator in where clause rendering.
    """
    with pytest.raises(bq_exceptions.SQLSyntaxParseException):
        mako_helpers.where(filters={"a__like": "x%"})


def test_add_cell_info_required_columns_none() -> None:
    """
    Return only required columns when input is None.
    """
    assert mako_helpers.add_cell_info_required_columns(columns=None) == ["c.id", "c.ingest_id"]


def test_add_cell_info_required_columns_adds_missing() -> None:
    """
    Prepend required columns when they are missing from the input list.
    """
    cols = ["c.foo", "c.bar"]
    result = mako_helpers.add_cell_info_required_columns(columns=cols)
    assert result[:2] == ["c.id", "c.ingest_id"]
    assert result[2:] == cols


def test_add_cell_info_required_columns_no_duplicates() -> None:
    """
    Preserve original list when it already contains both required columns.
    """
    cols = ["c.id", "c.ingest_id", "c.foo"]
    result = mako_helpers.add_cell_info_required_columns(columns=cols)
    assert result == cols


def test_remove_leading_alias() -> None:
    """
    Strip table aliases from fully qualified column names.
    """
    assert mako_helpers.remove_leading_alias(column_names=["t1.col1", "alias2.col2", "col3"]) == [
        "col1",
        "col2",
        "col3",
    ]


def test_select_basic_and_duplicates_preserved_case() -> None:
    """
    Render a select clause deduplicating columns while preserving case.
    """
    # Note: current implementation does not lowercase identifiers
    stmt = mako_helpers.select(column_names=["c.CELL", "i.DATASET", "c.CELL"])  # duplicates removed
    assert stmt == "select c.CELL, i.DATASET"


def test_select_empty_returns_star() -> None:
    """
    Return an asterisk when no columns are provided to select().
    """
    stmt = mako_helpers.select(column_names=[])
    assert stmt == "*"


def test_select_invalid_column_raises() -> None:
    """
    Raise an error when a column name contains more than one period.
    """
    with pytest.raises(ValueError):
        mako_helpers.select(column_names=["a.b.c"])  # too many periods
