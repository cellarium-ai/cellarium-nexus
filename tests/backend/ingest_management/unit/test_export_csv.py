from __future__ import annotations

import csv
from io import StringIO

import pytest

from cellarium.nexus.backend.cell_management import models as cm_models
from cellarium.nexus.backend.ingest_management import models as im_models
from cellarium.nexus.backend.ingest_management.utils.export_csv import export_model_queryset_to_csv


def _parse_response_csv(resp) -> list[list[str]]:
    content = resp.content.decode("utf-8")
    reader = csv.reader(StringIO(content))
    return [row for row in reader]


def test_export_model_queryset_to_csv_raises_on_empty_queryset() -> None:
    qs = cm_models.BigQueryDataset.objects.none()
    with pytest.raises(ValueError):
        export_model_queryset_to_csv(queryset=qs)


def test_export_model_queryset_to_csv_default_filename_and_header(default_dataset: cm_models.BigQueryDataset) -> None:
    qs = cm_models.BigQueryDataset.objects.filter(pk=default_dataset.pk)

    resp = export_model_queryset_to_csv(queryset=qs)

    # Filename header
    disp = resp["Content-Disposition"]
    assert disp.startswith("attachment; filename=bigquerydataset_export_")
    assert disp.endswith(".csv")

    # Header row should match model field names in order
    expected_fields = [f.name for f in cm_models.BigQueryDataset._meta.fields]
    rows = _parse_response_csv(resp)
    assert rows[0] == expected_fields

    # Data row contains the instance values
    # We just assert name/link show up, and id matches
    assert rows[1][expected_fields.index("id")] == str(default_dataset.pk)
    assert rows[1][expected_fields.index("name")] == default_dataset.name


def test_export_model_queryset_to_csv_exclude_fields(default_dataset: cm_models.BigQueryDataset) -> None:
    qs = cm_models.BigQueryDataset.objects.filter(pk=default_dataset.pk)

    resp = export_model_queryset_to_csv(queryset=qs, exclude_fields=["description"])
    rows = _parse_response_csv(resp)

    # Header should not include excluded field
    assert "description" not in rows[0]


def test_export_model_queryset_to_csv_foreign_key_id(default_dataset: cm_models.BigQueryDataset) -> None:
    ingest = im_models.IngestInfo.objects.create(bigquery_dataset=default_dataset)

    qs = im_models.IngestInfo.objects.filter(pk=ingest.pk)
    resp = export_model_queryset_to_csv(queryset=qs)
    rows = _parse_response_csv(resp)

    header = rows[0]
    data = rows[1]

    # Find FK column index
    idx = header.index("bigquery_dataset")
    assert data[idx] == str(default_dataset.pk)
