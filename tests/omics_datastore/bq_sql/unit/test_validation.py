import pytest

from cellarium.nexus.omics_datastore.bq_ops.bq_sql.validation import exceptions as val_exceptions
from cellarium.nexus.omics_datastore.bq_ops.bq_sql.validation import template_data_validator as v


def test_validate_not_empty_raises() -> None:
    """
    Raise a TemplateDataValidationError when value is empty.
    """
    with pytest.raises(val_exceptions.TemplateDataValidationError):
        v.validate_not_empty(value="")


def test_validate_sql_filter_eq_and_in() -> None:
    """
    Accept eq with primitives and in with homogeneous lists; reject heterogeneous lists.
    """
    # eq accepts primitives
    v.validate_sql_filter(filter_name="a__eq", filter_value=1)
    v.validate_sql_filter(filter_name="a__eq", filter_value=True)
    v.validate_sql_filter(filter_name="a__eq", filter_value=1.5)
    v.validate_sql_filter(filter_name="a__eq", filter_value="x")

    # in accepts homogeneous lists
    v.validate_sql_filter(filter_name="b__in", filter_value=[1, 2])
    v.validate_sql_filter(filter_name="c__in", filter_value=["a", "b"])

    with pytest.raises(val_exceptions.TemplateDataValidationError):
        v.validate_sql_filter(filter_name="c__in", filter_value=[1, "b"])  # heterogenous


def test_validate_sql_filter_unsupported_type_message() -> None:
    """
    Raise a TemplateDataValidationError for unsupported filter operator and include guidance in message.
    """
    with pytest.raises(val_exceptions.TemplateDataValidationError) as ei:
        v.validate_sql_filter(filter_name="a__like", filter_value="x%")
    # Message references eq/in availability currently
    assert "eq" in str(ei.value) and "in" in str(ei.value)


def test_validate_column_name_rules() -> None:
    """
    Validate column names accept one or zero periods and reject more than one or empty names.
    """
    v.validate_column_name(column_name="t.col")
    v.validate_column_name(column_name="col")
    with pytest.raises(ValueError):
        v.validate_column_name(column_name="a.b.c")
    with pytest.raises(ValueError):
        v.validate_column_name(column_name="")
