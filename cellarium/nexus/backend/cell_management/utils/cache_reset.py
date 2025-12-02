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
from django.conf import settings
from google.cloud import bigquery

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.utils import bigquery_utils
from cellarium.nexus.omics_datastore.bq_ops.data_operator import BigQueryDataOperator

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
    # Retrieve datasets
    datasets_qs = cell_models.OmicsDataset.objects.all().order_by("name")
    dataset_names: list[str] = [ds.name for ds in datasets_qs]

    bq_client = bigquery.Client(project=settings.GCP_PROJECT_ID)

    for ds in dataset_names:
        if not ds:
            continue

        # Create manager for this dataset
        operator = BigQueryDataOperator(
            client=bq_client,
            project=settings.GCP_PROJECT_ID,
            dataset=ds,
        )
        manager = bigquery_utils.OmicsCachedDataManager(
            operator=operator,
            cache_namespace=f"bq|{settings.GCP_PROJECT_ID}|{ds}",
        )

        # Warm count for empty filters
        try:
            manager.get_cached_count(filters_dict={})
            repopulated_items += 1
        except Exception as exc:  # noqa: BLE001 - continue warming others
            logger.warning(f"Failed to warm count cache for dataset '{ds}': {exc}")

        # Warm categorical columns
        try:
            categorical = manager.get_cached_categorical_obs_columns()
            repopulated_items += 1
        except Exception as exc:  # noqa: BLE001 - continue warming others
            logger.warning(f"Failed to warm categorical columns for dataset '{ds}': {exc}")
            categorical = set()

        # Warm distinct values for each categorical column
        for col in sorted(categorical):
            try:
                manager.get_cached_distinct_obs_values(column_name=col)
                repopulated_items += 1
            except Exception as exc:  # noqa: BLE001 - continue warming others
                logger.warning(f"Failed to warm distinct values for dataset '{ds}', column '{col}': {exc}")

    logger.info(f"Repopulated {repopulated_items} cache entries across {len(dataset_names)} datasets")
    return {"deleted": int(deleted_rows), "repopulated": int(repopulated_items)}
