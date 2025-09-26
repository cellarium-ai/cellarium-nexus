"""
Utility functions for resetting Django cache entries used by the admin UI.

Delegate the cache reset to cell-management utilities that understand the
BigQuery-backed filter caching and perform a targeted reset.
"""

import logging

from cellarium.nexus.backend.cell_management.utils import cache_reset as cm_cache_reset

logger = logging.getLogger(__name__)


def reset_cache_and_repopulate() -> list[str]:
    """
    Reset and repopulate cache for BigQuery-backed filters.

    Delegate to ``cell_management`` utility that handles targeted deletion
    of cache entries (e.g., ``countcache:*`` and ``bqcache:*``) and warms
    common entries to improve first-load performance.

    :raise RuntimeError: If reset fails at the database/cache layer
    :raise Exception: If BigQuery operations fail during repopulation

    :return: Summary of cache namespaces that were reset
    """
    logger.info("Resetting BigQuery-backed filter cache via cell_management utility")
    stats = cm_cache_reset.reset_cellinfo_filter_cache(repopulate=True)
    logger.info(f"Cache reset details: {stats}")
    # Return a simple summary of the namespaces we reset
    return ["countcache:*", "bqcache:*"]
