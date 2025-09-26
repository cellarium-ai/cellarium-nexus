from __future__ import annotations

import pytest

from cellarium.nexus.backend.cell_management import models


@pytest.mark.usefixtures("backend_cache_cleaner")
def test_feature_schema_roundtrip(default_dataset: models.BigQueryDataset) -> None:
    """
    Persist ``FeatureInfo`` objects and ensure schema relations resolve.

    :param default_dataset: Default dataset fixture
    """

    feature = models.FeatureInfo.objects.create(ensemble_id="ENSG000001", symbol="TP53")
    schema = models.FeatureSchema.objects.create(name="default-schema")
    schema.features.add(feature)

    retrieved = models.FeatureSchema.objects.get(name="default-schema")

    assert retrieved.features.count() == 1
    assert retrieved.features.first().ensemble_id == "ENSG000001"


def test_get_default_dataset_only_when_single(default_dataset: models.BigQueryDataset) -> None:
    """
    Evaluate ``BigQueryDatasetQuerySet.get_default_dataset`` behavior.

    :param default_dataset: Default dataset fixture
    """

    singleton = models.BigQueryDataset.objects.get_default_dataset()
    assert singleton == default_dataset

    models.BigQueryDataset.objects.create(name="secondary", description="secondary")

    assert models.BigQueryDataset.objects.get_default_dataset() is None
