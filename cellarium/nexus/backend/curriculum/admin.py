from typing import Any

from django.contrib import admin
from django.db import models
from django.http import HttpRequest
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter

from cellarium.nexus.backend.cell_management.admin.forms import CustomJSONEditorWidget
from cellarium.nexus.backend.curriculum.models import Curriculum
from cellarium.nexus.backend.shared.admin import GCSDisplayMixin


@admin.register(Curriculum)
class CurriculumAdmin(GCSDisplayMixin, ModelAdmin):
    """
    Admin interface for managing curriculums.

    Provides functionality to manage curriculum metadata and extraction settings.
    """

    list_display: list[str] = [
        "id",
        "name",
        "creator",
        "cell_count",
        "extract_bin_size",
        "status",
        "created_at",
    ]
    list_filter: list[str | tuple[str, type]] = [
        ("id", RangeNumericFilter),
        ("cell_count", RangeNumericFilter),
        ("extract_bin_size", RangeNumericFilter),
        ("creator", RelatedDropdownFilter),
        "status",
        "created_at",
    ]
    search_fields: list[str] = [
        "name",
        "creator__username",
    ]
    readonly_fields: list[str] = ["created_at", "status"]
    ordering: list[str] = ["-created_at"]

    list_gcs_viewable: list[dict[str, Any]] = [
        {"name": "extract_files_path", "is_directory": True},
        {"name": "metadata_file_path", "is_directory": False},
    ]

    def get_queryset(self, request: HttpRequest) -> models.QuerySet[Curriculum]:
        """
        Get the queryset for the admin view.

        Add select_related for better performance.

        :param request: The HTTP request

        :return: Queryset with related fields
        """
        return super().get_queryset(request).select_related("creator")

    def formfield_for_dbfield(self, db_field: models.Field, request: HttpRequest, **kwargs) -> Any:
        """
        Override formfield to use custom widgets for specific fields.

        Handle FilePathField and JSON fields with appropriate widgets.

        :param db_field: Database field
        :param request: HTTP request

        :return: Form field with appropriate widget
        """
        # Handle JSON editor widget
        if db_field.name == "filters_json":
            kwargs["widget"] = CustomJSONEditorWidget(
                attrs={
                    "height": "400px",
                    "width": "100%",
                },
                options={
                    "modes": ["tree", "code"],
                    "mode": "tree",
                    "search": True,
                    "sortObjectKeys": True,
                    "enableSort": False,
                    "enableTransform": False,
                },
            )
        return super().formfield_for_dbfield(db_field, request, **kwargs)
