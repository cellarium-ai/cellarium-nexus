from __future__ import annotations

import json

import pytest
from django.test import Client
from rest_framework.reverse import reverse

from cellarium.nexus.backend.cell_management import models as cm_models


@pytest.mark.usefixtures("backend_cache_cleaner")
@pytest.mark.django_db
def test_cellinfo_filters_count_invalid_json(client: Client, admin_user: object) -> None:
    client.force_login(user=admin_user)
    url = reverse("admin:cell_management_cellinfo_filters_count")

    resp = client.post(path=url, data=b"not-json", content_type="application/json")
    assert resp.status_code == 400
    payload = resp.json()
    assert payload.get("error") == "invalid_json"


@pytest.mark.usefixtures("backend_cache_cleaner")
@pytest.mark.django_db
def test_cellinfo_filters_count_missing_dataset(client: Client, admin_user: object) -> None:
    client.force_login(user=admin_user)
    url = reverse("admin:cell_management_cellinfo_filters_count")

    resp = client.post(path=url, data=json.dumps({"dataset": "", "filters": {}}), content_type="application/json")
    assert resp.status_code == 400
    payload = resp.json()
    assert payload.get("error") == "invalid_dataset"


@pytest.mark.usefixtures("backend_cache_cleaner", "bigquery_cached_manager_stub")
@pytest.mark.django_db
def test_cellinfo_filters_count_success(
    client: Client, admin_user: object, default_dataset: cm_models.BigQueryDataset
) -> None:
    client.force_login(user=admin_user)
    url = reverse("admin:cell_management_cellinfo_filters_count")

    body = {"dataset": default_dataset.name, "filters": {}}
    resp = client.post(path=url, data=json.dumps(body), content_type="application/json")
    assert resp.status_code == 200
    assert resp.json()["count"] == 42


@pytest.mark.usefixtures("backend_cache_cleaner", "bigquery_cached_manager_stub")
@pytest.mark.django_db
def test_cellinfo_admin_view_renders(
    client: Client, admin_user: object, default_dataset: cm_models.BigQueryDataset
) -> None:
    client.force_login(user=admin_user)
    url = reverse("admin:cell_management_cellinfo_changelist")

    resp = client.get(path=url)
    assert resp.status_code == 200
