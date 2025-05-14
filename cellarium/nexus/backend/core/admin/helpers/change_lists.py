import logging

from django.contrib.admin.views.main import ChangeList
from django.core.paginator import Paginator
from django.db.models import QuerySet
from django.http import HttpRequest

from cellarium.nexus.backend.core.admin.utils.caching_utils import get_cached_count

logger = logging.getLogger(__name__)


class CachingPaginator(Paginator):
    """
    A custom paginator that uses cached count operations.

    This paginator overrides the count property to use cached counts
    instead of direct database queries.
    """

    @property
    def count(self) -> int:
        """
        Get the count of objects, using caching.

        This property overrides the default count behavior to use caching.

        :return: The count of objects
        """
        if not hasattr(self, "_count"):
            self._count = get_cached_count(self.object_list)
        return self._count


class CachingChangeList(ChangeList):
    """
    A ChangeList subclass that caches count operations to improve performance.

    This class overrides the get_results method to use cached counts instead of
    performing expensive count queries on each request.
    """

    def get_queryset(self, request: HttpRequest) -> QuerySet:
        """
        Get the queryset for the changelist.

        :param request: The HTTP request object

        :return: The filtered queryset
        """
        return super().get_queryset(request)

    def get_results(self, request: HttpRequest) -> None:
        """
        Get the results for the changelist page with cached counts.

        This method overrides the parent implementation to use cached counts.

        :param request: The HTTP request object

        :raise: Exception: Any exceptions that might occur during queryset operations

        :return: None
        """
        self.queryset = self.get_queryset(request)

        paginator = CachingPaginator(self.queryset, self.list_per_page)
        self.paginator = paginator

        if not self.show_all:
            try:
                self.result_list = paginator.get_page(self.page_num).object_list
            except Exception:
                self.result_list = paginator.get_page(1).object_list
        else:
            self.result_list = self.queryset

        self.result_count = get_cached_count(self.queryset)

        if self.is_popup and self.list_select_related:
            full_result_queryset = self.root_queryset
        else:
            full_result_queryset = self.queryset

        self.full_result_count = get_cached_count(full_result_queryset)

        self.can_show_all = self.full_result_count <= self.list_max_show_all
        self.multi_page = self.result_count > self.list_per_page
