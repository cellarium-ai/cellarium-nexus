import hashlib
import logging

from django.conf import settings
from django.core.cache import cache
from django.db.models import QuerySet

logger = logging.getLogger(__name__)


def get_cached_count(queryset: QuerySet, timeout: int = settings.COUNT_CACHE_TTL_SECONDS) -> int:
    """
    Get count of queryset with caching.

    Creates a cache key based on the SQL query and parameters, then checks if
    the count is already cached. If not, performs the count operation and caches
    the result.

    :param queryset: The queryset to count
    :param timeout: Cache timeout in seconds, defaults to 30 days

    :raise: ValueError: If queryset is not a valid QuerySet instance

    :return: Count of objects in the queryset
    """
    sql, params = queryset.query.sql_with_params()
    raw = sql + str(params)
    key = f"countcache_{hashlib.md5(raw.encode()).hexdigest()}"

    # Try to get the cached count
    cached = cache.get(key)
    if cached is not None:
        return cached

    count = queryset.count()

    cache.set(key, count, timeout)
    return count
