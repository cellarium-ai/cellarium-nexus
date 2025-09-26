from django.contrib.auth import get_user_model

from cellarium.nexus.backend import curriculum
from cellarium.nexus.backend.curriculum.api import serializers


def test_curriculum_serializer_roundtrip(db) -> None:
    """
    Assert CurriculumSerializer maps creator_id <-> creator and respects read-only fields.
    """
    User = get_user_model()
    user = User.objects.create_user(username="alice", password="x")

    # Create via serializer using creator_id
    payload = {
        "name": "my-curriculum",
        "creator_id": user.id,
        "cell_count": 100,
        "extract_bin_size": 500,
        "extract_files_path": "gs://bucket/extracts/",
        "metadata_file_path": "gs://bucket/meta.json",
        "filters_json": {"organism__eq": "Homo sapiens"},
        "status": curriculum.models.Curriculum.STATUS_PREPARE,
    }
    ser = serializers.CurriculumSerializer(data=payload)
    assert ser.is_valid(), ser.errors
    instance = ser.save()

    assert instance.creator_id == user.id
    assert instance.name == "my-curriculum"

    # Read-only fields should be present but not writable: mutate and ensure ignored on update
    original_id = instance.id
    original_created_at = instance.created_at

    readback = serializers.CurriculumSerializer(instance)
    data = readback.data
    assert data["id"] == original_id
    assert data["creator_id"] == user.id

    # Attempt to update read-only fields; serializer should ignore them
    update_ser = serializers.CurriculumSerializer(
        instance,
        data={
            "name": "my-curriculum-updated",
            "creator_id": user.id,
            "id": "should-be-ignored",
            "created_at": "2000-01-01T00:00:00Z",
        },
        partial=True,
    )
    assert update_ser.is_valid(), update_ser.errors
    updated = update_ser.save()

    assert updated.id == original_id
    assert updated.created_at == original_created_at
    assert updated.name == "my-curriculum-updated"
