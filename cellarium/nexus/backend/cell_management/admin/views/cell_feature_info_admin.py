"""
Admin module for cell feature information management.
"""

from django.contrib import admin
from django.http import HttpRequest
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter

from cellarium.nexus.backend.cell_management.admin.filters import TagDropdownFilter
from cellarium.nexus.backend.cell_management.models import CellFeatureInfo


@admin.register(CellFeatureInfo)
class CellFeatureInfoAdmin(ModelAdmin):
    """
    Admin interface for managing cell feature information.

    Provides functionality to filter and search cell feature data, including numeric range filters
    for the id field. Prevents direct creation and editing of cell feature info instances.
    """

    list_display = (
        "id",
        "ensemble_id",
        "symbol",
        "biotype",
        "is_filtered",
        "reference",
        "ingest",
        "tag",
    )
    search_fields = (
        "id",
        "ensemble_id",
        "symbol",
        "biotype",
        "reference",
        "tag",
    )
    list_filter = (
        ("id", RangeNumericFilter),
        "is_filtered",
        TagDropdownFilter,
        ("ingest", RelatedDropdownFilter),
    )
    ordering = ("-id",)
    readonly_fields = ("id",)

    def has_add_permission(self, request: HttpRequest) -> bool:
        """
        Disable the ability to create new cell feature info instances directly.

        :param request: The HTTP request

        :return: False to prevent direct creation
        """
        return False

    def has_change_permission(self, request: HttpRequest, obj=None) -> bool:
        """
        Disable the ability to edit cell feature info instances.

        :param request: The HTTP request
        :param obj: The object being changed

        :return: False to prevent editing
        """
        return False
