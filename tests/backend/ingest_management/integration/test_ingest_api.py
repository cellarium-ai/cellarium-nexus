from __future__ import annotations

import json

import pytest
from django.core.cache import cache
from django.test import Client
from rest_framework.reverse import reverse

from cellarium.nexus.backend.cell_management import models as cm_models
from cellarium.nexus.backend.ingest_management import models as im_models


@pytest.mark.django_db
def test_ingest_create_retrieve_and_update(client: Client, default_dataset: cm_models.BigQueryDataset) -> None:
    """
    Create an ingest via API, retrieve it, then update fields via PATCH.
    """
    # Create
    create_url = reverse("ingest-create")
    payload = {"bigquery_dataset": default_dataset.name}
    resp = client.post(path=create_url, data=json.dumps(payload), content_type="application/json")
    assert resp.status_code in (200, 201)
    created = resp.json()
    ingest_id = created["id"]

    # Retrieve
    detail_url = reverse("ingest-retrieve-update", kwargs={"id": ingest_id})
    resp = client.get(path=detail_url)
    assert resp.status_code == 200
    data = resp.json()
    assert data["bigquery_dataset"] == default_dataset.name

    # Update metadata_extra and status via PATCH
    patch_payload = {"metadata_extra": {"foo": "bar"}, "status": im_models.IngestInfo.STATUS_SUCCEEDED}
    resp = client.patch(path=detail_url, data=json.dumps(patch_payload), content_type="application/json")
    assert resp.status_code == 200
    updated = resp.json()
    assert updated["metadata_extra"] == {"foo": "bar"}
    assert updated["status"] == im_models.IngestInfo.STATUS_SUCCEEDED


@pytest.mark.django_db
def test_reserve_indexes_cell_info_db_side_effects(client: Client, default_dataset: cm_models.BigQueryDataset) -> None:
    """
    Reserve indexes for cell_info and assert IndexTracking row is created/updated.
    """

    url = reverse("reserve-indexes-cell-info")
    payload = {"bigquery_dataset": default_dataset.name, "batch_size": 100}

    # First reservation should create tracking row and allocate 1..100
    resp = client.post(path=url, data=json.dumps(payload), content_type="application/json")
    assert resp.status_code == 200
    assert resp.json() == {"index_start": 1, "index_end": 100}
    tracking = im_models.IndexTracking.objects.get(bigquery_dataset=default_dataset, resource_key="cell_info")
    assert tracking.largest_index == 100

    # Second reservation should continue from 101
    payload = {"bigquery_dataset": default_dataset.name, "batch_size": 50}
    resp = client.post(path=url, data=json.dumps(payload), content_type="application/json")
    assert resp.status_code == 200
    assert resp.json() == {"index_start": 101, "index_end": 150}
    tracking.refresh_from_db()
    assert tracking.largest_index == 150


@pytest.mark.django_db
def test_reserve_indexes_feature_info_db_side_effects(
    client: Client, default_dataset: cm_models.BigQueryDataset
) -> None:
    """
    Reserve indexes for feature_info and assert IndexTracking row is created/updated.
    """
    url = reverse("reserve-indexes-feature-info")
    payload = {"bigquery_dataset": default_dataset.name, "batch_size": 50}

    resp = client.post(path=url, data=json.dumps(payload), content_type="application/json")
    assert resp.status_code == 200
    assert resp.json() == {"index_start": 1, "index_end": 50}
    tracking = im_models.IndexTracking.objects.get(bigquery_dataset=default_dataset, resource_key="feature_info")
    assert tracking.largest_index == 50


@pytest.mark.django_db
def test_validation_report_item_create(client: Client, admin_user: object) -> None:
    """
    Create a ValidationReportItem via API for an existing report.
    """
    client.force_login(user=admin_user)

    report = im_models.ValidationReport.objects.create(creator=admin_user)
    url = reverse("validation-report-item-create")
    payload = {
        "report_id": report.id,
        "input_file_gcs_path": "gs://bucket/path/file.h5ad",
        "validator_name": "validate_schema",
        "is_valid": True,
        "message": "ok",
    }
    resp = client.post(path=url, data=json.dumps(payload), content_type="application/json")
    assert resp.status_code in (200, 201)

    # Ensure persisted
    assert im_models.ValidationReportItem.objects.filter(report=report, validator_name="validate_schema").exists()


@pytest.mark.django_db
@pytest.mark.usefixtures("bigquery_cached_manager_stub")
def test_reset_cache_endpoint_e2e(client: Client, admin_user: object) -> None:
    """
    POST reset-cache and assert cache side effects and response envelope.
    """
    client.force_login(user=admin_user)

    # Seed namespaced cache keys and ensure present
    cache.clear()
    cache.set("countcache:seed1", 1)
    cache.set("bqcache:seed2", ["a", "b"])
    cache.set("other:untouched", 1)

    url = reverse("reset-cache")
    resp = client.post(path=url)
    assert resp.status_code == 200
    out = resp.json()
    assert out["status"] == "success"
    # Function returns summary patterns; order not guaranteed
    assert set(out["repopulated_keys"]) == {"countcache:*", "bqcache:*"}

    # Ensure targeted keys are deleted; unrelated keys remain
    assert cache.get("countcache:seed1") is None
    assert cache.get("bqcache:seed2") is None
    assert cache.get("other:untouched") == 1
