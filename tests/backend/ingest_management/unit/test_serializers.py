"""
Unit tests for ingest management API serializers.
"""

import pytest
from rest_framework.exceptions import ValidationError

from cellarium.nexus.backend.cell_management import models as cell_management_models
from cellarium.nexus.backend.ingest_management import models as ingest_models
from cellarium.nexus.backend.ingest_management.api import serializers


@pytest.fixture()
def ingest(default_dataset: cell_management_models.BigQueryDataset) -> ingest_models.IngestInfo:
    """
    Create ingest info bound to the provided default dataset.

    :param default_dataset: Default BigQuery dataset fixture
    :return: Newly created ingest info instance
    """
    return ingest_models.IngestInfo.objects.create(bigquery_dataset=default_dataset)


def test_ingest_from_avro_validate_happy_path(ingest: ingest_models.IngestInfo) -> None:
    """
    Validate serializer on happy path and ensure ingest relation is resolved.
    """
    s = serializers.IngestFromAvroSerializer(data={"stage_dir": "/tmp/x", "ingest_id": ingest.id})
    assert s.is_valid(), s.errors
    data = s.validated_data
    assert data["ingest"].id == ingest.id


def test_ingest_from_avro_validate_status_check(ingest: ingest_models.IngestInfo) -> None:
    """
    Raise when ingest is already in a terminal succeeded state.

    :raise: ValidationError
    """
    # Mark as succeeded to trigger error
    ingest.status = ingest_models.IngestInfo.STATUS_SUCCEEDED
    ingest.save(update_fields=["status"])

    s = serializers.IngestFromAvroSerializer(data={"stage_dir": "/tmp/x", "ingest_id": ingest.id})
    with pytest.raises(ValidationError):
        s.is_valid(raise_exception=True)


def test_ingest_from_avro_nonexistent_id() -> None:
    """
    Raise when ingest id does not correspond to any existing record.

    :raise: ValidationError
    """
    s = serializers.IngestFromAvroSerializer(data={"stage_dir": "/tmp/x", "ingest_id": 999999})
    with pytest.raises(ValidationError):
        s.is_valid(raise_exception=True)


def test_reserve_indexes_serializer_slug_field(default_dataset: cell_management_models.BigQueryDataset) -> None:
    """
    Validate slug field resolves dataset by name and returns model instance.
    """
    s = serializers.ReserveIndexesSerializer(data={"bigquery_dataset": default_dataset.name, "batch_size": 10})
    assert s.is_valid(), s.errors
    assert s.validated_data["bigquery_dataset"].id == default_dataset.id
