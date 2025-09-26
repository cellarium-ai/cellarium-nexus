"""
Unit tests for ingest management API serializers.
"""

from typing import Any, Callable

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


def _set_status(ingest: ingest_models.IngestInfo, *, status: str) -> None:
    """
    Persist ingest status update for parametrized tests.

    :param ingest: Ingest info instance to update
    :param status: Status value to apply
    """
    ingest.status = status
    ingest.save(update_fields=["status"])


@pytest.mark.parametrize(
    "payload_factory, setup_case",
    [
        pytest.param(
            lambda ingest: {"stage_dir": "/tmp/x", "ingest_id": ingest.id},
            lambda ingest: _set_status(ingest, status=ingest_models.IngestInfo.STATUS_SUCCEEDED),
            id="status-succeeded",
        ),
        pytest.param(
            lambda ingest: {"stage_dir": "/tmp/x", "ingest_id": 999999},
            lambda ingest: None,
            id="missing-ingest",
        ),
    ],
)
def test_ingest_from_avro_validation_errors(
    ingest: ingest_models.IngestInfo,
    payload_factory: Callable[[ingest_models.IngestInfo], dict[str, Any]],
    setup_case: Callable[[ingest_models.IngestInfo], None] | None,
) -> None:
    """
    Validate serializer error scenarios for ingest lookup failures.

    :param ingest: Default ingest info instance provided by fixture
    :param payload_factory: Callable that builds serializer payload for the scenario
    :param setup_case: Callable that mutates ingest state prior to serializer validation

    :raise: ValidationError
    """
    if setup_case is not None:
        setup_case(ingest)

    payload = payload_factory(ingest)
    serializer = serializers.IngestFromAvroSerializer(data=payload)

    with pytest.raises(ValidationError):
        serializer.is_valid(raise_exception=True)


def test_reserve_indexes_serializer_slug_field(default_dataset: cell_management_models.BigQueryDataset) -> None:
    """
    Validate slug field resolves dataset by name and returns model instance.
    """
    s = serializers.ReserveIndexesSerializer(data={"bigquery_dataset": default_dataset.name, "batch_size": 10})
    assert s.is_valid(), s.errors
    assert s.validated_data["bigquery_dataset"].id == default_dataset.id
