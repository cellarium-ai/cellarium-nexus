from django.contrib import admin
from django.core.cache import cache
from django.db import models
from django.http import HttpRequest
from unfold.contrib.filters.admin import DropdownFilter, MultipleDropdownFilter

CACHE_KEY_FORMAT = "generic_dropdown_filter:{parameter_name}"


class CachedLookupsMixin:
    """
    Mixin for caching lookup values in dropdown filters.

    This mixin provides caching functionality for dropdown filters,
    allowing them to store and retrieve lookup values from the cache.
    """

    field_name: str = None

    def _get_cache_key(self) -> str:
        """
        Get the cache key for this filter.

        :raise: NotImplementedError if not implemented by subclass

        :return: Cache key string
        """
        raise NotImplementedError("Subclasses must implement this method.")

    @classmethod
    def _get_unique_values_from_db_or_cache(cls, cache_key: str, model: models.Model) -> list[str]:
        """
        Get unique values from cache or database.

        Retrieve values from cache if available, otherwise query the database
        and store the results in cache.

        :param cache_key: The cache key to use
        :param model: The model to query if cache is empty

        :return: List of unique values
        """
        values = cache.get(cache_key)

        if values is None:
            values = model.objects.values_list(cls.field_name, flat=True).distinct()
            cache.set(cache_key, values, timeout=None)

        return values

    def lookups(self, request: HttpRequest, model_admin: admin.ModelAdmin) -> list[tuple[str, str]]:
        """
        Get lookup values for the filter.

        :param request: The HTTP request
        :param model_admin: The model admin instance

        :return: List of tuples containing lookup values
        """
        cache_key = self._get_cache_key()
        values = self._get_unique_values_from_db_or_cache(cache_key=cache_key, model=model_admin.model)
        return [(v, v) for v in values if v not in ["", None]]


class GenericDropdownFilter(CachedLookupsMixin, DropdownFilter):
    """Generic dropdown filter with caching support.

    This filter provides a single-select dropdown with cached options.
    """

    title = None
    parameter_name = None
    field_name = None

    def _get_cache_key(self) -> str:
        """
        Get the cache key for this filter.

        :return: Cache key string
        """
        return CACHE_KEY_FORMAT.format(parameter_name=self.parameter_name)

    def queryset(self, request: HttpRequest, queryset: models.QuerySet) -> models.QuerySet:
        """
        Filter the queryset based on the selected value.

        :param request: The HTTP request
        :param queryset: The queryset to filter

        :return: Filtered queryset
        """
        if self.value() not in ["", None]:
            return queryset.filter(**{self.field_name: self.value()})
        return queryset


class GenericMultiDropdownFilter(CachedLookupsMixin, MultipleDropdownFilter):
    """Generic multi-select dropdown filter with caching support.

    This filter provides a multi-select dropdown with cached options.
    """

    title = None
    parameter_name = None
    field_name = None

    def _get_cache_key(self) -> str:
        """
        Get the cache key for this filter.

        :return: Cache key string
        """
        return CACHE_KEY_FORMAT.format(parameter_name=self.parameter_name)

    def queryset(self, request: HttpRequest, queryset: models.QuerySet) -> models.QuerySet:
        """
        Filter the queryset based on the selected values.

        :param request: The HTTP request
        :param queryset: The queryset to filter

        :return: Filtered queryset
        """
        values = self.value()
        if values:
            return queryset.filter(**{f"{self.field_name}__in": values})
        return queryset
