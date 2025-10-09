from __future__ import annotations

import json

import pytest

from cellarium.nexus.backend.cell_management import models as cm_models
from cellarium.nexus.backend.cell_management.admin import forms as admin_forms


def _make_feature_schema() -> cm_models.FeatureSchema:
    fs = cm_models.FeatureSchema.objects.create(name="fs-test")
    cm_models.FeatureInfo.objects.create(
        ensemble_id="E1",
        symbol="S1",
        feature_schema=fs,
    )
    return fs


def _valid_form_data(*, dataset: cm_models.BigQueryDataset) -> dict:
    return {
        "feature_schema": _make_feature_schema().pk,
        "name": "uniq-name",
        "extract_bin_size": 1000,
        "bigquery_dataset": dataset.pk,
        "obs_columns": ["organism"],
        "categorical_column_count_limit": 2000,
        # read-only widget, but still posted/cleaned
        "filters": json.dumps({}),
    }


@pytest.mark.usefixtures("dummy_gcs_client")
def test_extract_curriculum_form_clean_name_db_duplicate(
    default_dataset: cm_models.BigQueryDataset, admin_user
) -> None:
    # Create a curriculum with the same name in DB
    from cellarium.nexus.backend.curriculum import models as cur_models

    cur_models.Curriculum.objects.create(
        name="dup", creator=admin_user, cell_count=1, extract_bin_size=1, extract_bin_count=1
    )

    data = _valid_form_data(dataset=default_dataset)
    data["name"] = "dup"

    form = admin_forms.ExtractCurriculumForm(data=data)
    assert not form.is_valid()
    assert "name" in form.errors


@pytest.mark.usefixtures("dummy_gcs_client")
def test_extract_curriculum_form_clean_name_gcs_duplicate(
    default_dataset: cm_models.BigQueryDataset, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Pretend GCS contains a curriculum prefix with that name
    monkeypatch.setattr(
        "cellarium.nexus.backend.cell_management.admin.forms.check_curriculum_exists",
        lambda name: True,
        raising=True,
    )

    data = _valid_form_data(dataset=default_dataset)
    data["name"] = "exists-in-bucket"

    form = admin_forms.ExtractCurriculumForm(data=data)
    assert not form.is_valid()
    assert "name" in form.errors


@pytest.mark.usefixtures("dummy_gcs_client")
def test_extract_curriculum_form_clean_filters_parsing(default_dataset: cm_models.BigQueryDataset) -> None:
    data = _valid_form_data(dataset=default_dataset)
    data["filters"] = json.dumps({"organism__eq": ["human"]})

    form = admin_forms.ExtractCurriculumForm(data=data)
    assert form.is_valid()
    assert form.cleaned_data["filters"] == {"organism__eq": ["human"]}


def test_custom_json_editor_widget_format_value() -> None:
    w = admin_forms.CustomJSONEditorWidget()
    out = w.format_value('{"a":1}')
    assert out == {"a": 1}
