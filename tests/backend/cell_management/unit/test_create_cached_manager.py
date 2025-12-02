"""
Unit tests for _create_cached_manager_for_dataset factory function.
"""

from __future__ import annotations

import pytest

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.admin.views.cell_info_admin import (
    _create_cached_manager_for_dataset,
)
from cellarium.nexus.backend.cell_management.utils import bigquery_utils
from cellarium.nexus.omics_datastore.bq_ops.data_operator import BigQueryDataOperator
from cellarium.nexus.omics_datastore.soma_ops.data_operator import TileDBSOMADataOperator


@pytest.mark.django_db
def test_create_cached_manager_bigquery_backend(default_dataset: cell_models.OmicsDataset) -> None:
    """
    Assert BigQuery backend creates BigQueryDataOperator.

    :param default_dataset: Default BigQuery dataset fixture
    """
    # Ensure dataset is BigQuery backend
    default_dataset.backend = cell_models.OmicsDatasetBackend.BIGQUERY
    default_dataset.save()

    manager = _create_cached_manager_for_dataset(dataset_name=default_dataset.name)

    assert isinstance(manager, bigquery_utils.OmicsCachedDataManager)
    assert isinstance(manager.operator, BigQueryDataOperator)
    assert "bq|" in manager.cache_namespace


@pytest.mark.django_db
def test_create_cached_manager_soma_backend() -> None:
    """
    Assert TileDB SOMA backend creates TileDBSOMADataOperator.
    """
    # Create a SOMA dataset
    soma_dataset = cell_models.OmicsDataset.objects.create(
        name="test-soma-dataset",
        backend=cell_models.OmicsDatasetBackend.TILEDB_SOMA,
        uri="gs://test-bucket/soma-experiment",
    )

    manager = _create_cached_manager_for_dataset(dataset_name=soma_dataset.name)

    assert isinstance(manager, bigquery_utils.OmicsCachedDataManager)
    assert isinstance(manager.operator, TileDBSOMADataOperator)
    assert manager.operator.experiment_uri == "gs://test-bucket/soma-experiment"
    assert "soma|" in manager.cache_namespace


@pytest.mark.django_db
def test_create_cached_manager_soma_without_uri_raises() -> None:
    """
    Assert TileDB SOMA backend without URI raises ValueError.
    """
    # Create a SOMA dataset without URI
    soma_dataset = cell_models.OmicsDataset.objects.create(
        name="test-soma-no-uri",
        backend=cell_models.OmicsDatasetBackend.TILEDB_SOMA,
        uri=None,
    )

    with pytest.raises(ValueError, match="has no URI configured"):
        _create_cached_manager_for_dataset(dataset_name=soma_dataset.name)


@pytest.mark.django_db
def test_create_cached_manager_nonexistent_dataset_raises() -> None:
    """
    Assert non-existent dataset raises DoesNotExist.
    """
    with pytest.raises(cell_models.OmicsDataset.DoesNotExist):
        _create_cached_manager_for_dataset(dataset_name="nonexistent-dataset")
