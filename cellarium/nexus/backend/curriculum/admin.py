from django.contrib import admin
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter

from cellarium.nexus.backend.curriculum.models import Curriculum


@admin.register(Curriculum)
class CurriculumAdmin(ModelAdmin):
    """
    Admin interface for managing curriculums.

    Provides functionality to manage curriculum metadata and extraction settings.
    """

    list_display = [
        "id",
        "creator",
        "cell_count",
        "extract_bin_size",
        "extract_files_dir",
        "metadata_file_path",
        "created_at",
    ]
    list_filter = [
        ("id", RangeNumericFilter),
        ("cell_count", RangeNumericFilter),
        ("extract_bin_size", RangeNumericFilter),
        ("creator", RelatedDropdownFilter),
        "created_at",
    ]
    search_fields = [
        "extract_files_dir",
        "metadata_file_path",
        "creator__username",
    ]
    readonly_fields = ["created_at"]
    ordering = ["-created_at"]
    fieldsets = (
        (
            None,
            {
                "fields": (
                    "creator",
                    "cell_count",
                    "extract_bin_size",
                    "extract_files_dir",
                    "metadata_file_path",
                    "filters_json",
                    "created_at",
                )
            },
        ),
    )

    def get_queryset(self, request):
        """
        Get the queryset for the admin view.

        Add select_related for better performance.

        :param request: The HTTP request
        :return: Queryset with related fields
        """
        return super().get_queryset(request).select_related("creator")
