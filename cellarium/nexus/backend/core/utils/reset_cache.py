"""
Utility functions for resetting and repopulating Django cache.

This module provides functionality to clear all cache values and 
repopulate the cache with cached filters from the admin interface.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

from django.apps import apps
from django.core.cache import cache

from cellarium.nexus.backend.cell_management.admin.filters import (
    CellTypeDropdownFilter,
    DiseaseDropdownFilter,
    DonorDropdownFilter,
    OrganismDropdownFilter,
    SexDropdownFilter,
    SuspensionTypeDropdownFilter,
    TagDropdownFilter,
)
from cellarium.nexus.backend.core.admin.filters.cache_dropdown_filters import CACHE_KEY_FORMAT

logger = logging.getLogger(__name__)


def reset_cache_and_repopulate() -> list[str]:
    """
    Clear all cache values and repopulate cache with cached filters.

    This function performs a complete cache reset and then repopulates
    the cache with all the cached filters used in the admin interface.
    It requires a working database connection as it will query the database
    to populate the filter caches.

    :raise: RuntimeError if any filter repopulation fails

    :return: List of cache keys that were repopulated
    """
    logger.info("Clearing all cache values")
    cache.clear()

    filter_classes = [
        SexDropdownFilter,
        DonorDropdownFilter,
        SuspensionTypeDropdownFilter,
        OrganismDropdownFilter,
        TagDropdownFilter,
        CellTypeDropdownFilter,
        DiseaseDropdownFilter,
    ]

    try:
        cell_info_model = apps.get_model("cell_management", "CellInfo")
    except LookupError as e:
        logger.error(f"Failed to get CellInfo model: {e}")
        raise RuntimeError(f"Failed to get CellInfo model: {e}")

    repopulated_keys = []

    def repopulate_filter_cache(filter_class: Any) -> tuple[str, int]:
        cache_key = CACHE_KEY_FORMAT.format(parameter_name=filter_class.parameter_name)
        values = filter_class._get_unique_values_from_db_or_cache(cache_key=cache_key, model=cell_info_model)

        return cache_key, len(values)

    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_filter = {
            executor.submit(repopulate_filter_cache, filter_class): filter_class for filter_class in filter_classes
        }
        for future in as_completed(future_to_filter):
            filter_class = future_to_filter[future]
            try:
                cache_key, values_count = future.result()
                logger.info(f"Repopulated cache for {filter_class.__name__} with {values_count} values")
                repopulated_keys.append(cache_key)
            except Exception as e:
                logger.error(f"Failed to repopulate cache for {filter_class.__name__}: {e}")
                raise RuntimeError(f"Failed to repopulate cache for {filter_class.__name__}: {e}")

    logger.info(f"Successfully repopulated {len(repopulated_keys)} cache keys")
    return repopulated_keys
