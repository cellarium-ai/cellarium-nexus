import hashlib
import json
import logging
from typing import Any

from django.conf import settings
from django.core.cache import cache
from google.cloud import bigquery

from cellarium.nexus.omics_datastore import bq_ops
from cellarium.nexus.omics_datastore.bq_ops import constants as bq_constants

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


class BigQueryCachedDataManager:
    """
    Provide cached access to BigQuery cell counts using a shared client.

    :param project_id: Optional GCP project ID; if not provided, use settings.GCP_PROJECT_ID

    :return: None
    """

    def __init__(self, project_id: str = settings.GCP_PROJECT_ID) -> None:
        self.project_id = project_id
        self.bq_client = bigquery.Client(project=self.project_id)

    def get_cached_count_bq(
        self,
        *,
        dataset_name: str,
        filters_dict: dict[str, Any],
        timeout: int = settings.COUNT_CACHE_TTL_SECONDS,
    ) -> int:
        """
        Get cell count from BigQuery and cache the result using a stable key.

        The cache key derives from project ID, dataset name, and a canonical JSON
        serialization of filters to ensure stability regardless of dict order.

        :param dataset_name: Name of the BigQuery dataset to query
        :param filters_dict: Dictionary of filter statements for the count query
        :param timeout: Cache TTL in seconds

        :raise: google.api_core.exceptions.GoogleAPICallError

        :return: Integer count of cells returned by the filter
        """

        serialized_filters = _canonical_json(filters_dict)
        key_seed = f"v1|{self.project_id}|{dataset_name}|{serialized_filters}"
        key = _hash_to_key(prefix="countcache", seed=key_seed)

        cached = cache.get(key)
        if cached is not None:
            return int(cached)

        operator = bq_ops.BigQueryDataOperator(
            client=self.bq_client,
            project=self.project_id,
            dataset=dataset_name,
        )
        count = operator.count_cells(filter_statements=filters_dict)
        cache.set(key, int(count), timeout)
        return int(count)

    def get_cached_categorical_columns_bq(
        self,
        *,
        dataset_name: str,
        table_name: str = bq_constants.BQ_CELL_INFO_TABLE_NAME,
        distinct_threshold: int = settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
        timeout: int = settings.SUGGEST_CACHE_TTL_SECONDS,
    ) -> set[str]:
        """
        Get categorical string columns for a base table and cache the result.

        :param dataset_name: BigQuery dataset name
        :param table_name: Base table name to analyze
        :param distinct_threshold: Threshold for distinct count to classify as categorical
        :param timeout: Optional cache TTL in seconds

        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

        :return: Set of categorical column names
        """
        threshold = distinct_threshold

        key_seed = f"v1|categorical_columns|{self.project_id}|{dataset_name}|{table_name}|{threshold}"
        key = _hash_to_key(prefix="bqcache", seed=key_seed)

        cached = cache.get(key)
        if cached is not None:
            return set(cached)

        operator = bq_ops.BigQueryDataOperator(
            client=self.bq_client,
            project=self.project_id,
            dataset=dataset_name,
        )
        categorical = operator.get_categorical_string_columns(
            table_name=table_name,
            threshold=int(threshold),
        )
        cache.set(key, list(categorical), timeout)
        return categorical

    def get_cached_distinct_values_bq(
        self,
        *,
        dataset_name: str,
        column_name: str,
        table_name: str = bq_constants.BQ_CELL_INFO_TABLE_NAME,
        limit: int = settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
        timeout: int = settings.SUGGEST_CACHE_TTL_SECONDS,
    ) -> list[str]:
        """
        Get distinct values for a base table column and cache the result.

        :param dataset_name: BigQuery dataset name
        :param column_name: Column to fetch distinct values for
        :param table_name: Base table name to query
        :param limit: Maximum number of values to return
        :param timeout: Cache TTL in seconds

        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

        :return: List of distinct values as strings
        """
        key_seed = f"v1|distinct_values|{self.project_id}|{dataset_name}|{table_name}|{column_name}|{limit}"
        key = _hash_to_key(prefix="bqcache", seed=key_seed)

        cached = cache.get(key)
        if cached is not None:
            return list(cached)

        operator = bq_ops.BigQueryDataOperator(
            client=self.bq_client,
            project=self.project_id,
            dataset=dataset_name,
        )
        values = operator.get_distinct_values(
            table_name=table_name, column_name=column_name, limit=limit, dataset=dataset_name
        )
        cache.set(key, list(values), timeout)
        return values
