import hashlib
import logging
from typing import Any

from django.conf import settings
from django.core.cache import cache
from google.cloud import bigquery

from cellarium.nexus.omics_datastore.bq_ops import BigQueryDataOperator

logger = logging.getLogger(__name__)


class BigQueryCachedDataManager:
    def __init__(self, bigquery_dataset_name: str):
        bq_client = bigquery.Client(project=settings.GCP_PROJECT_ID)
        self.bq_data_operator = BigQueryDataOperator(
            client=bq_client, project=settings.GCP_PROJECT_ID, dataset=bigquery_dataset_name
        )

    def get_cached_count_bq(self, filters_dict: dict[str, Any], timeout: int = settings.COUNT_CACHE_TTL_SECONDS) -> int:
        """
        Get count of queryset with caching.

        Creates a cache key based on the SQL query and parameters, then checks if
        the count is already cached. If not, performs the count operation and caches
        the result.

        :param filters_dict: Dictionary of filters to apply to the queryset
        :param bigquery_dataset_name: Name of the BigQuery dataset
        :param timeout: Cache timeout in seconds, defaults to 30 days

        :raise: ValueError: If queryset is not a valid QuerySet instance

        :return: Count of objects in the queryset
        """
        raw = str(filters_dict)
        key = f"countcache_{hashlib.md5(raw.encode()).hexdigest()}"

        # Try to get the cached count
        cached = cache.get(key)
        if cached is not None:
            return cached

        count = self.bq_data_operator.count_cells(filter_statements=filters_dict)

        cache.set(key, count, timeout)
        return count
