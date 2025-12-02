from __future__ import annotations

from typing import Any

import pytest
from django.test import Client

from cellarium.nexus.backend.application import app_storages_backend


def test_root_redirects_to_admin(client: Client) -> None:
    """
    Assert root URL redirects to admin dashboard.

    :param client: Django test client fixture
    """

    response = client.get(path="/")

    assert response.status_code == 302
    assert response.headers["Location"].endswith("/admin/")


@pytest.mark.usefixtures("default_dataset", "omics_cached_manager_stub")
def test_custom_admin_urls_accessible(
    admin_client: Client,
) -> None:
    """
    Ensure custom admin URLs are registered and render successfully.

    :param admin_client: Django admin client fixture
    """

    response = admin_client.get(path="/admin/cell_management/cellinfo/")

    assert response.status_code == 200
    assert "datasets" in response.context


@pytest.mark.parametrize(
    "input_name,expected",
    [
        pytest.param("folder/file.txt", "gs://test-bucket/folder/file.txt", id="simple"),
        pytest.param("/leading/slash.txt", "gs://test-bucket/leading/slash.txt", id="leading-slash"),
    ],
)
@pytest.mark.usefixtures("dummy_gcs_client")
def test_custom_google_cloud_storage_path(
    settings: Any,
    input_name: str,
    expected: str,
) -> None:
    """
    Compute canonical ``gs://`` path using custom storage backend.

    :param settings: Django settings fixture
    :param input_name: File name passed to storage backend
    :param expected: Expected gs path output
    """

    settings.GS_BUCKET_NAME = "test-bucket"
    storage_backend = app_storages_backend.CustomGoogleCloudStorage()

    result = storage_backend.gs_path(name=input_name)

    assert result == expected
