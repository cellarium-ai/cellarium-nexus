import hashlib
import json
import logging
from typing import Any

from django.conf import settings
from django.core.cache import cache

from cellarium.nexus.omics_datastore.protocols import DataOperatorProtocol

logger = logging.getLogger(__name__)


def _canonical_json(data: dict[str, Any]) -> str:
    """
    Produce deterministic JSON for use in cache keys.

    :param data: Mapping to serialize

    :return: Canonical JSON string
    """
    return json.dumps(data, sort_keys=True, separators=(",", ":"))


def _hash_to_key(prefix: str, seed: str) -> str:
    """
    Build a cache key by hashing a seed with MD5 and prefixing a namespace.

    :param prefix: Namespace prefix for the key (e.g., "countcache", "bqcache")
    :param seed: Seed string to hash

    :return: Namespaced cache key
    """
    return f"{prefix}:{hashlib.md5(seed.encode('utf-8')).hexdigest()}"


class OmicsCachedDataManager:
    """
    Provide cached access to omics data using a data operator.

    This class wraps any operator conforming to DataOperatorProtocol and adds
    caching for expensive operations like cell counts and distinct value queries.

    :param operator: Data operator implementing DataOperatorProtocol
    :param cache_namespace: Unique namespace for cache keys (e.g., dataset name)
    """

    def __init__(
        self,
        *,
        operator: DataOperatorProtocol,
        cache_namespace: str,
    ) -> None:
        self.operator = operator
        self.cache_namespace = cache_namespace

    def get_cached_count(
        self,
        *,
        filters_dict: dict[str, Any],
        timeout: int = settings.COUNT_CACHE_TTL_SECONDS,
    ) -> int:
        """
        Get cell count and cache the result using a stable key.

        The cache key derives from namespace and a canonical JSON serialization
        of filters to ensure stability regardless of dict order.

        :param filters_dict: Dictionary of filter statements for the count query
        :param timeout: Cache TTL in seconds

        :raise Exception: If operator query fails

        :return: Integer count of cells returned by the filter
        """
        serialized_filters = _canonical_json(filters_dict)
        key_seed = f"v1|{serialized_filters}"
        key = _hash_to_key(prefix=f"countcache|{self.cache_namespace}", seed=key_seed)

        cached = cache.get(key)
        if cached is not None:
            return int(cached)

        count = self.operator.count_cells(filter_statements=filters_dict)
        cache.set(key, int(count), timeout)
        return int(count)

    def get_cached_categorical_obs_columns(
        self,
        *,
        distinct_threshold: int = settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
        timeout: int = settings.SUGGEST_CACHE_TTL_SECONDS,
    ) -> set[str]:
        """
        Get categorical string columns for obs/cell_info and cache the result.

        :param distinct_threshold: Threshold for distinct count to classify as categorical
        :param timeout: Cache TTL in seconds

        :raise Exception: If operator query fails

        :return: Set of categorical column names
        """
        key_seed = f"v1|{distinct_threshold}"
        key = _hash_to_key(
            prefix=f"omicscache|{self.cache_namespace}|categorical_obs_columns",
            seed=key_seed,
        )

        cached = cache.get(key)
        if cached is not None:
            return set(cached)

        categorical = self.operator.get_categorical_obs_columns(threshold=distinct_threshold)
        cache.set(key, list(categorical), timeout)
        return categorical

    def get_cached_distinct_obs_values(
        self,
        *,
        column_name: str,
        timeout: int = settings.SUGGEST_CACHE_TTL_SECONDS,
    ) -> list[str]:
        """
        Get distinct values for an obs/cell_info column and cache the result.

        :param column_name: Column to fetch distinct values for
        :param timeout: Cache TTL in seconds

        :raise Exception: If operator query fails

        :return: List of distinct values as strings
        """
        key_seed = f"v1|{column_name}"
        key = _hash_to_key(
            prefix=f"omicscache|{self.cache_namespace}|distinct_obs_values",
            seed=key_seed,
        )

        cached = cache.get(key)
        if cached is not None:
            return list(cached)

        values = self.operator.get_distinct_obs_values(column_name=column_name)
        cache.set(key, list(values), timeout)
        return values
