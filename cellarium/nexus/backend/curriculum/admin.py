from django.contrib import admin
from django.utils.html import format_html
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter

from cellarium.nexus.backend.cell_management.admin.forms import CustomJSONEditorWidget
from cellarium.nexus.backend.curriculum.models import Curriculum


@admin.register(Curriculum)
class CurriculumAdmin(ModelAdmin):
    """
    Admin interface for managing curriculums.

    Provides functionality to manage curriculum metadata and extraction settings.
    """

    list_display = [
        "id",
        "name",
        "creator",
        "cell_count",
        "extract_bin_size",
        "status",
        "display_extract_files_dir",
        "display_metadata_file_path",
        "created_at",
    ]
    list_filter = [
        ("id", RangeNumericFilter),
        ("cell_count", RangeNumericFilter),
        ("extract_bin_size", RangeNumericFilter),
        ("creator", RelatedDropdownFilter),
        "status",
        "created_at",
    ]
    search_fields = [
        "name",
        "extract_files_dir",
        "metadata_file_path",
        "creator__username",
    ]
    readonly_fields = ["created_at", "status", "display_extract_files_dir", "display_metadata_file_path"]
    ordering = ["-created_at"]
    fieldsets = (
        (
            None,
            {
                "fields": (
                    "name",
                    "creator",
                    "cell_count",
                    "extract_bin_size",
                    "status",
                    "display_extract_files_dir",
                    "display_metadata_file_path",
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
        :raise: None
        """
        return super().get_queryset(request).select_related("creator")

    def formfield_for_dbfield(self, db_field, request, **kwargs):
        """
        Override formfield to use custom JSON widget for filters_json.

        :param db_field: Database field
        :param request: HTTP request
        :return: Form field with appropriate widget
        :raise: None
        """
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

    def display_extract_files_dir(self, obj):
        """
        Display extract files directory in a more readable format.

        :param obj: Curriculum instance
        :return: HTML formatted path display
        :raise: None
        """
        if not obj.extract_files_dir:
            return "-"

        path_parts = obj.extract_files_dir.split("/")
        if len(path_parts) <= 2:
            return obj.extract_files_dir

        return format_html(
            '<span title="{}"><i class="fas fa-folder"></i> {}/…/{}</span>',
            obj.extract_files_dir,
            path_parts[0],
            path_parts[-1],
        )

    display_extract_files_dir.short_description = "Extract Files Directory"
    display_extract_files_dir.admin_order_field = "extract_files_dir"

    def display_metadata_file_path(self, obj):
        """
        Display metadata file path in a more readable format.

        :param obj: Curriculum instance
        :return: HTML formatted path display
        :raise: None
        """
        if not obj.metadata_file_path:
            return "-"

        path_parts = obj.metadata_file_path.split("/")
        filename = path_parts[-1]
        directory = "/".join(path_parts[:-1])

        return format_html(
            '<span title="{}"><i class="fas fa-file"></i> {}</span>',
            obj.metadata_file_path,
            filename,
        )

    display_metadata_file_path.short_description = "Metadata File"
    display_metadata_file_path.admin_order_field = "metadata_file_path"
