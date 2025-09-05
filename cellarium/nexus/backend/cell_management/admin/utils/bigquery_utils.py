import hashlib
import json
import logging
from typing import Any

from django.conf import settings
from django.core.cache import cache
from google.cloud import bigquery

from cellarium.nexus.omics_datastore import bq_ops

logger = logging.getLogger(__name__)


class BigQueryCachedDataManager:
    """
    Provide cached access to BigQuery cell counts using a shared client.

    :param project_id: Optional GCP project ID; if not provided, use settings.GCP_PROJECT_ID

    :return: None
    """

    def __init__(self, project_id: str | None = None) -> None:
        self.project_id = project_id or settings.GCP_PROJECT_ID
        self.bq_client = bigquery.Client(project=self.project_id)

    def get_cached_count_bq(
        self,
        *,
        dataset_name: str,
        filters_dict: dict[str, Any],
        timeout: int | None = None,
    ) -> int:
        """
        Get cell count from BigQuery and cache the result using a stable key.

        The cache key derives from project ID, dataset name, and a canonical JSON
        serialization of filters to ensure stability regardless of dict order.

        :param dataset_name: Name of the BigQuery dataset to query
        :param filters_dict: Dictionary of filter statements for the count query
        :param timeout: Optional cache TTL in seconds; if not provided, use 
            settings.COUNT_CACHE_TTL_SECONDS (defaulting to 30 days if missing)

        :raise: google.api_core.exceptions.GoogleAPICallError

        :return: Integer count of cells returned by the filter
        """
        # Resolve timeout at call time to avoid binding settings at import time
        default_ttl = 60 * 60 * 24 * 30  # 30 days
        ttl = timeout if timeout is not None else getattr(settings, "COUNT_CACHE_TTL_SECONDS", default_ttl)

        # Build a stable cache key
        serialized_filters = json.dumps(filters_dict, sort_keys=True, separators=(",", ":"))
        key_seed = f"v1|{self.project_id}|{dataset_name}|{serialized_filters}"
        key = f"countcache:{hashlib.md5(key_seed.encode('utf-8')).hexdigest()}"

        # Try cache first
        cached = cache.get(key)
        if cached is not None:
            return int(cached)

        operator = bq_ops.BigQueryDataOperator(
            client=self.bq_client,
            project=self.project_id,
            dataset=dataset_name,
        )
        count = operator.count_cells(filter_statements=filters_dict)
        cache.set(key, int(count), ttl)
        return int(count)
