from __future__ import annotations

import pytest
from django.core.cache import cache

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.utils.cache_reset import (
    reset_cellinfo_filter_cache,
)


def test_reset_cellinfo_filter_cache_delete_only() -> None:
    """
    Seed DatabaseCache with namespaced keys and assert targeted deletion.

    :return: None
    """
    # Ensure clean state (backend_cache_cleaner is autouse but keep explicit here)
    cache.clear()

    # Seed keys in both targeted namespaces and a control key
    cache.set("countcache:seed1", 1)
    cache.set("bqcache:seed2", ["a", "b"])
    cache.set("other:untouched", 1)

    result = reset_cellinfo_filter_cache(repopulate=False)

    assert result["repopulated"] == 0
    # Expect both targeted keys to be deleted
    assert result["deleted"] == 2

    assert cache.get("countcache:seed1") is None
    assert cache.get("bqcache:seed2") is None
    # Non-targeted keys should remain
    assert cache.get("other:untouched") == 1


@pytest.mark.usefixtures("omics_cached_manager_stub")
@pytest.mark.django_db
def test_reset_cellinfo_filter_cache_repopulate_basic(default_dataset: models.OmicsDataset) -> None:
    """
    Repopulate using the cached manager stub; expect 2 items per dataset:
    1 for counts; 1 for categorical columns; no distincts when categories are empty.

    :param default_dataset: Default BigQuery dataset fixture

    :return: None
    """
    # No pre-seeded keys
    cache.clear()

    result = reset_cellinfo_filter_cache(repopulate=True)

    dataset_count = models.OmicsDataset.objects.count()
    assert result["deleted"] == 0
    assert result["repopulated"] == 2 * dataset_count


@pytest.mark.usefixtures("omics_cached_manager_stub")
@pytest.mark.django_db
def test_reset_cellinfo_filter_cache_repopulate_with_distincts(
    default_dataset: models.OmicsDataset, monkeypatch: pytest.MonkeyPatch
) -> None:
    """
    Repopulate when categorical columns exist; expect 3 items per dataset:
    1 for counts; 1 for categorical columns; 1 for distinct values for that column.

    :param default_dataset: Default BigQuery dataset fixture
    :param monkeypatch: Pytest monkeypatch fixture

    :return: None
    """
    # Import the stub class installed by the fixture and adjust return values
    from tests.backend.fixtures.mocks import OmicsCachedManagerStub

    monkeypatch.setattr(
        OmicsCachedManagerStub,
        "get_cached_categorical_obs_columns",
        lambda self, **kwargs: {"organism"},
        raising=True,
    )
    monkeypatch.setattr(
        OmicsCachedManagerStub,
        "get_cached_distinct_obs_values",
        lambda self, **kwargs: ["Homo sapiens"],
        raising=True,
    )

    cache.clear()

    result = reset_cellinfo_filter_cache(repopulate=True)

    dataset_count = models.OmicsDataset.objects.count()
    assert result["deleted"] == 0
    assert result["repopulated"] == 3 * dataset_count
