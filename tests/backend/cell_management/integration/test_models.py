from __future__ import annotations

import pytest

from cellarium.nexus.backend.cell_management import models


@pytest.mark.usefixtures("backend_cache_cleaner")
def test_feature_schema_roundtrip(default_dataset: models.BigQueryDataset) -> None:
    """
    Persist ``FeatureInfo`` objects and ensure schema relations resolve.

    :param default_dataset: Default dataset fixture
    """

    schema = models.FeatureSchema.objects.create(name="default-schema")
    _ = models.FeatureInfo.objects.create(
        ensemble_id="ENSG000001",
        symbol="TP53",
        feature_schema=schema,
    )

    retrieved = models.FeatureSchema.objects.get(name="default-schema")

    assert retrieved.features.count() == 1
    assert retrieved.features.first().ensemble_id == "ENSG000001"


@pytest.mark.usefixtures("backend_cache_cleaner")
def test_feature_schema_cascade_delete(default_dataset: models.BigQueryDataset) -> None:
    """
    Verify that deleting a schema cascades to delete its features.

    :param default_dataset: Default dataset fixture
    """

    schema = models.FeatureSchema.objects.create(name="test-schema")
    feature1 = models.FeatureInfo.objects.create(
        ensemble_id="ENSG000001",
        symbol="TP53",
        feature_schema=schema,
    )
    feature2 = models.FeatureInfo.objects.create(
        ensemble_id="ENSG000002",
        symbol="EGFR",
        feature_schema=schema,
    )

    assert models.FeatureInfo.objects.count() == 2

    # Delete the schema
    schema.delete()

    # Verify features are also deleted
    assert models.FeatureInfo.objects.count() == 0
    assert not models.FeatureInfo.objects.filter(pk=feature1.pk).exists()
    assert not models.FeatureInfo.objects.filter(pk=feature2.pk).exists()


def test_get_default_dataset_only_when_single(default_dataset: models.BigQueryDataset) -> None:
    """
    Evaluate ``BigQueryDatasetQuerySet.get_default_dataset`` behavior.

    :param default_dataset: Default dataset fixture
    """

    singleton = models.BigQueryDataset.objects.get_default_dataset()
    assert singleton == default_dataset

    models.BigQueryDataset.objects.create(name="secondary", description="secondary")

    assert models.BigQueryDataset.objects.get_default_dataset() is None
