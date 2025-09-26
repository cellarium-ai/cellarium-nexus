"""
Utilities to reset and repopulate cache entries for BigQuery-backed filters.

This module targets cache keys produced by BigQuery-backed filter utilities:
- countcache:* (cell counts per dataset and filters)
- bqcache:* (categorical columns and distinct values)
"""

from __future__ import annotations

import logging

import django.conf as dj_conf
import django.db as dj_db

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.utils import bigquery_utils
from cellarium.nexus.omics_datastore.bq_ops import constants as bq_constants

logger = logging.getLogger(__name__)


def reset_cellinfo_filter_cache(*, repopulate: bool = True) -> dict[str, int]:
    """
    Reset cache entries for BigQuery-backed filter data and optionally repopulate them.

    This function targets cache keys with the following prefixes:
    - ``countcache:`` for cached counts
    - ``bqcache:`` for cached categorical columns and distinct values

    Assumes the default cache backend is DatabaseCache and performs a targeted
    SQL deletion of rows from the cache table.

    When ``repopulate`` is ``True``, this function warms the cache by computing:
    - per-dataset cell counts with empty filters
    - per-dataset categorical string columns for the base cell info table
    - per-dataset distinct values for each categorical column (limited by settings)

    :param repopulate: Whether to repopulate caches after reset

    :raise RuntimeError: If database operations fail during reset
    :raise Exception: If BigQuery operations fail during repopulation

    :return: Mapping with counts of deleted rows and repopulated items
    """
    deleted_rows = 0
    repopulated_items = 0

    # Targeted deletion for DatabaseCache backend by SQL
    caches = dj_conf.settings.CACHES or {}
    default_cache = caches.get("default", {})
    location = default_cache.get("LOCATION")
    if not location:
        raise RuntimeError("DatabaseCache LOCATION is not configured for default cache.")

    with dj_db.connection.cursor() as cursor:
        # Pull keys for logging/count (best effort) and delete them
        select_sql = f"SELECT cache_key FROM {location} WHERE cache_key LIKE %s OR cache_key LIKE %s"
        delete_sql = f"DELETE FROM {location} WHERE cache_key LIKE %s OR cache_key LIKE %s"
        # Use wildcard patterns to match DatabaseCache versioned keys (e.g., ":1:countcache:...")
        cursor.execute(select_sql, ["%countcache:%", "%bqcache:%"])
        keys = [row[0] for row in cursor.fetchall()] or []
        deleted_rows = len(keys)

        cursor.execute(delete_sql, ["%countcache:%", "%bqcache:%"])

    logger.info(f"Deleted {deleted_rows} cache rows with prefixes 'countcache:' and 'bqcache:' from table '{location}'")

    if not repopulate:
        return {"deleted": int(deleted_rows), "repopulated": int(repopulated_items)}

    # Repopulate caches by warming commonly used entries
    manager = bigquery_utils.BigQueryCachedDataManager()

    # Retrieve datasets
    datasets_qs = cell_models.BigQueryDataset.objects.all().order_by("name")
    dataset_names: list[str] = [ds.name for ds in datasets_qs]

    for ds in dataset_names:
        if not ds:
            continue
        # Warm count for empty filters
        try:
            manager.get_cached_count_bq(dataset_name=ds, filters_dict={})
            repopulated_items += 1
        except Exception as exc:  # noqa: BLE001 - continue warming others
            logger.warning(f"Failed to warm count cache for dataset '{ds}': {exc}")

        # Warm categorical columns
        try:
            categorical = manager.get_cached_categorical_columns_bq(
                dataset_name=ds,
                table_name=bq_constants.BQ_CELL_INFO_TABLE_NAME,
            )
            repopulated_items += 1
        except Exception as exc:  # noqa: BLE001 - continue warming others
            logger.warning(f"Failed to warm categorical columns for dataset '{ds}': {exc}")
            categorical = set()

        # Warm distinct values for each categorical column
        for col in sorted(categorical):
            try:
                manager.get_cached_distinct_values_bq(
                    dataset_name=ds,
                    column_name=col,
                    table_name=bq_constants.BQ_CELL_INFO_TABLE_NAME,
                )
                repopulated_items += 1
            except Exception as exc:  # noqa: BLE001 - continue warming others
                logger.warning(f"Failed to warm distinct values for dataset '{ds}', column '{col}': {exc}")

    logger.info(f"Repopulated {repopulated_items} cache entries across {len(dataset_names)} datasets")
    return {"deleted": int(deleted_rows), "repopulated": int(repopulated_items)}
