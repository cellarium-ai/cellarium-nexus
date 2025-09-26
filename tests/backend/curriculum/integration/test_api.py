from __future__ import annotations

import json

from django.test import Client
from rest_framework.reverse import reverse

from cellarium.nexus.backend.curriculum import models


def test_curriculum_create_and_retrieve(
    client: Client,
    curriculum_factory,
    admin_user,
) -> None:
    """
    Exercise Curriculum API create and retrieve endpoints.

    :param client: Django client fixture
    :param curriculum_factory: Factory creating curriculum instances
    """

    payload = {
        "name": "integration-curriculum",
        "creator_id": admin_user.id,
        "cell_count": 123,
        "filters_json": {"species": "Mus musculus"},
    }

    create_url = reverse(viewname="curriculum-list")
    client.force_login(user=admin_user)
    response = client.post(
        path=create_url,
        data=json.dumps(payload),
        content_type="application/json",
    )

    assert response.status_code == 201

    detail_url = reverse(viewname="curriculum-detail", kwargs={"name": "integration-curriculum"})
    detail_response = client.get(path=detail_url)

    assert detail_response.status_code == 200
    assert detail_response.json()["name"] == "integration-curriculum"


def test_curriculum_progress_update(client: Client, curriculum_factory) -> None:
    """
    Update curriculum status and ensure persisted state matches response.

    :param client: Django client fixture
    :param curriculum_factory: Factory creating curriculum instances
    """

    curriculum = curriculum_factory("progress-curriculum")
    client.force_login(user=curriculum.creator)

    detail_url = reverse(viewname="curriculum-detail", kwargs={"name": curriculum.name})
    payload = {"status": models.Curriculum.STATUS_SUCCEEDED}

    response = client.patch(path=detail_url, data=payload, content_type="application/json")

    assert response.status_code == 200

    curriculum.refresh_from_db()
    assert curriculum.status == models.Curriculum.STATUS_SUCCEEDED
