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

    pre_count = models.Curriculum.objects.count()

    payload = {
        "name": "integration-curriculum",
        "creator_id": admin_user.id,
        "cell_count": 123,
        "filters_json": {"species": "Mus musculus"},
    }

    create_url = reverse(viewname="curriculum-list")

    response = client.post(
        path=create_url,
        data=json.dumps(payload),
        content_type="application/json",
    )

    assert response.status_code == 201

    # Assert DB side effects: object created with expected fields
    assert models.Curriculum.objects.count() == pre_count + 1
    obj = models.Curriculum.objects.get(name="integration-curriculum")
    assert obj.creator_id == admin_user.id
    assert obj.cell_count == 123
    assert obj.filters_json == {"species": "Mus musculus"}
    assert obj.status == models.Curriculum.STATUS_PREPARE

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
    pre_count = models.Curriculum.objects.count()

    detail_url = reverse(viewname="curriculum-detail", kwargs={"name": curriculum.name})
    payload = {"status": models.Curriculum.STATUS_SUCCEEDED}

    response = client.patch(path=detail_url, data=payload, content_type="application/json")

    assert response.status_code == 200

    curriculum.refresh_from_db()
    assert curriculum.status == models.Curriculum.STATUS_SUCCEEDED
    # Count should remain unchanged on update
    assert models.Curriculum.objects.count() == pre_count
