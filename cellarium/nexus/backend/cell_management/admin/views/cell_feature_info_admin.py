"""
Admin module for cell feature information management.
"""

from django.contrib import admin
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import ChoicesDropdownFilter, RangeNumericFilter, RelatedDropdownFilter

from cellarium.nexus.backend.cell_management.admin.filters import MultiValueTextFilter
from cellarium.nexus.backend.cell_management.models import CellFeatureInfo


@admin.register(CellFeatureInfo)
class CellFeatureInfoAdmin(ModelAdmin):
    """
    Admin interface for managing cell feature information.

    Provides functionality to filter and search cell feature data, including numeric range filters
    for the id field.
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
        ("biotype", MultiValueTextFilter),
        ("is_filtered", ChoicesDropdownFilter),
        ("tag", MultiValueTextFilter),
        ("ensemble_id", MultiValueTextFilter),
        ("symbol", MultiValueTextFilter),
        ("ingest", RelatedDropdownFilter),
    )
    ordering = ("-id",)
    readonly_fields = ("id",)
