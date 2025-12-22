import pytest

from cellarium.nexus.omics_datastore.bq_ops.bq_sql import constants as constants_module
from cellarium.nexus.omics_datastore.bq_ops.bq_sql import template_data as td_module
from cellarium.nexus.omics_datastore.bq_ops.bq_sql.validation import exceptions as val_exceptions


def test_template_data_minimal() -> None:
    """
    Construct TemplateData with minimal arguments and verify stored keys.
    """
    tdata = td_module.TemplateData(project="proj")
    assert tdata.data[constants_module.TemplateDataDictNames.PROJECT] == "proj"
    assert constants_module.TemplateDataDictNames.DATASET not in tdata.data


def test_template_data_full_and_other_kwargs() -> None:
    """
    Construct TemplateData with dataset, select, filters, and other kwargs.
    """
    tdata = td_module.TemplateData(
        project="p",
        dataset="d",
        select=["c.id"],
        filters={"c.organism__eq": "Homo sapiens"},
        extra_flag=True,
        limit=100,
    )
    assert tdata.data[constants_module.TemplateDataDictNames.PROJECT] == "p"
    assert tdata.data[constants_module.TemplateDataDictNames.DATASET] == "d"
    assert tdata.data[constants_module.TemplateDataDictNames.SELECT] == ["c.id"]
    assert tdata.data[constants_module.TemplateDataDictNames.FILTERS]["c.organism__eq"] == "Homo sapiens"
    assert tdata.data["extra_flag"] is True
    assert tdata.data["limit"] == 100


def test_template_data_empty_project_raises() -> None:
    """
    Raise a TemplateDataValidationError when project is empty.
    """
    with pytest.raises(val_exceptions.TemplateDataValidationError):
        td_module.TemplateData(project="")


def test_template_data_invalid_filters_type_message() -> None:
    """
    Raise a TemplateDataValidationError for unsupported filter operator with informative message.
    """
    with pytest.raises(val_exceptions.TemplateDataValidationError) as ei:
        td_module.TemplateData(project="p", filters={"x__like": "abc%"})
    # Current error messaging mentions eq/in specifically
    assert "eq" in str(ei.value) and "in" in str(ei.value)


def test_template_data_filter_types_supported_passthrough() -> None:
    """
    Accept gt/gte/lt/lte filters and store them unchanged.
    """
    # Filters like gt/gte/lt/lte are allowed by validator (no extra checks performed)
    tdata = td_module.TemplateData(project="p", filters={"n__gt": 5, "m__lte": 3})
    assert tdata.data[constants_module.TemplateDataDictNames.FILTERS]["n__gt"] == 5
