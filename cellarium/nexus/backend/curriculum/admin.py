from django.conf import settings
from django.contrib import admin
from django.utils.html import format_html
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter

from cellarium.nexus.backend.cell_management.admin.forms import CustomJSONEditorWidget
from cellarium.nexus.backend.curriculum.models import Curriculum
from cellarium.nexus.shared.utils import gcp


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
        "display_extract_files_path",
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
        "creator__username",
    ]
    readonly_fields = ["created_at", "status", "extract_files_path", "metadata_file_path"]
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
                    "extract_files_path",
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
        :raise: None
        """
        return super().get_queryset(request).select_related("creator")

    def formfield_for_dbfield(self, db_field, request, **kwargs):
        """
        Override formfield to use custom widgets for specific fields.
        
        Handle FilePathField and JSON fields with appropriate widgets.

        :param db_field: Database field
        :param request: HTTP request
        
        :raise: None
        
        :return: Form field with appropriate widget
        """
        # Handle FilePathField to prevent errors when path doesn't exist
        if db_field.name == "extract_files_path":
            from django.forms import CharField
            from django.forms.widgets import TextInput
            return CharField(
                required=False,
                widget=TextInput(attrs={"readonly": "readonly"}),
                help_text=db_field.help_text,
                label=db_field.verbose_name,
                initial=db_field.value_from_object(db_field.model) if request.resolver_match.kwargs.get("object_id") else "",
            )
            
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

    # Use Django settings for GCS bucket names
    
    def _get_gcp_console_link(self, path: str, bucket_name: str = None, is_file: bool = False) -> tuple[str, str, str]:
        """
        Generate a GCP console link for a path, handling both full GCS paths and relative paths.
        
        :param path: Path to file or directory (can be full GCS path or relative path)
        :param bucket_name: Optional bucket name to use if path is relative
        :param is_file: Whether the path points to a file (True) or directory (False)
        
        :raise: ValueError - If the path is invalid and cannot be converted to a GCS path
        
        :return: Tuple containing (console_url, display_name, full_path)
        """
        if not path:
            return "", "", path
            
        # Handle full GCS paths (gs://bucket-name/path)
        if path.startswith("gs://"):
            try:
                bucket, file_path = gcp.get_bucket_name_and_file_path_from_gc_path(path)
                full_path = path
            except ValueError:
                return "", "", path
        # Handle relative paths
        else:
            # Use provided bucket name or get from Django settings
            bucket = bucket_name or settings.BUCKET_NAME_PRIVATE
            file_path = path.lstrip("/")
            full_path = f"gs://{bucket}/{file_path}"
        
        # Create a GCP console URL for the bucket/object
        # Use _details path only for files, not for directories
        if is_file:
            console_url = f"https://console.cloud.google.com/storage/browser/_details/{bucket}/{file_path}"
        else:
            console_url = f"https://console.cloud.google.com/storage/browser/{bucket}/{file_path}"
        
        # Create a display name (shortened version of the path)
        display_name = "View in Bucket"
        return console_url, display_name, full_path
    
    def display_extract_files_path(self, obj):
        """
        Display extract files path as a clickable link to GCP console.
        
        :param obj: Curriculum instance
        
        :raise: None
        
        :return: HTML formatted link to GCP console
        """
        if not obj.extract_files_path:
            return "-"
            
        # Get the path from the object
        path = str(obj.extract_files_path)
        
        # This is typically a directory, not a file
        # Generate the console link
        console_url, display_name, full_path = self._get_gcp_console_link(path, is_file=False)
        
        if not console_url:
            return path
            
        return format_html(
            '<a href="{}" target="_blank" title="{}"><i class="fas fa-folder"></i> {}</a>',
            console_url,
            full_path,
            display_name
        )
        
    display_extract_files_path.short_description = "Extract Files Path"
    display_extract_files_path.admin_order_field = "extract_files_path"

    def display_metadata_file_path(self, obj):
        """
        Display metadata file path as a clickable link to GCP console.
        
        :param obj: Curriculum instance
        
        :raise: None
        
        :return: HTML formatted link to GCP console
        """
        if not obj.metadata_file_path:
            return "-"
            
        # Get the path from the object
        file_path = str(obj.metadata_file_path)
        
        # For local files or non-GCS paths, just show the filename
        if file_path.startswith("/") or ("://" in file_path and not file_path.startswith("gs://")):
            path_parts = file_path.split("/")
            filename = path_parts[-1]
            return format_html(
                '<span title="{}"><i class="fas fa-file"></i> {}</span>',
                file_path,
                filename
            )
        
        # This is a file, so use the _details path
        # Generate the console link
        console_url, display_name, full_path = self._get_gcp_console_link(file_path, is_file=True)
        
        if not console_url:
            return file_path
            
        return format_html(
            '<a href="{}" target="_blank" title="{}"><i class="fas fa-file"></i> {}</a>',
            console_url,
            full_path,
            display_name
        )
        
    display_metadata_file_path.short_description = "Metadata File"
    display_metadata_file_path.admin_order_field = "metadata_file_path"
