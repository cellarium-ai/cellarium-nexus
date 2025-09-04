import logging
from typing import Any

from django.contrib.admin.views.main import ChangeList
from django.core.paginator import Paginator

from cellarium.nexus.backend.core.admin.utils.caching_utils import get_cached_count_bq

logger = logging.getLogger(__name__)


class BigQueryCountPaginator(Paginator):
    def __init__(
        self,
        object_list,
        per_page,
        orphans=0,
        allow_empty_first_page=True,
        bq_filters: dict[str, Any] | None = None,
        bigquery_dataset_name: str | None = None,
    ):
        super().__init__(
            object_list=object_list, per_page=per_page, orphans=orphans, allow_empty_first_page=allow_empty_first_page
        )
        self.bq_filters = bq_filters
        self.bigquery_dataset_name = bigquery_dataset_name

    @property
    def count(self) -> int:
        if self.bigquery_dataset_name is None:
            return 300

        return get_cached_count_bq(filters_dict=self.bq_filters, bigquery_dataset_name=self.bigquery_dataset_name)

    def page(self, number):
        # Return an empty page
        return super().page(1)
