"""
Admin mixins for cell management.
"""

from typing import Any

from django.conf import settings
from django.utils.html import format_html, mark_safe

from cellarium.nexus.shared.utils import gcp


class GCSDisplayMixin:
    """
    Add GCS file viewing capabilities to Django admin classes.

    This mixin adds functionality to display GCS file paths as clickable links in the Django admin
    interface, with appropriate icons for files and directories. It also adds view buttons to form
    fields for easy access to GCS resources.

    To use this mixin, subclass it in your ModelAdmin class and define a list_gcs_viewable attribute
    with the names of fields that should have GCS viewing capabilities.
    """

    list_gcs_viewable: list[str | dict[str, Any]] = []

    def __init_subclass__(cls, **kwargs) -> None:
        """
        Initialize the subclass by setting up GCS viewable fields.

        Process the list_gcs_viewable attribute and create display methods for each field,
        then add these methods to the list_display attribute.

        :raise: AttributeError: If the subclass has improperly configured list_gcs_viewable.
        """
        super().__init_subclass__(**kwargs)

        if not hasattr(cls, "list_display"):
            cls.list_display = []

        cls._gcs_view_fields_meta: dict[str, bool] = {}

        for item in cls.list_gcs_viewable:
            if isinstance(item, str):
                name = item
                is_directory = False
            else:
                name = item["name"]
                is_directory = item.get("is_directory", False)

            method_name = f"display_{name}"

            def make_display_func(field_name: str = name, is_dir: bool = is_directory) -> callable:
                """
                Create a display function for a GCS viewable field.

                :param field_name: The name of the field to display.
                :param is_dir: Whether the field represents a directory.

                :return: A function that renders the field as a clickable link.
                """

                def display_func(self, obj: Any) -> str:
                    """
                    Render a GCS field as a clickable link.

                    :param obj: The model instance being displayed.

                    :return: HTML markup for displaying the field.
                    """
                    value = getattr(obj, field_name)
                    if not value:
                        return "-"
                    value = str(value)

                    if value.startswith("/") or ("://" in value and not value.startswith("gs://")):
                        filename = value.split("/")[-1]
                        return format_html('<span title="{}"><i class="fas fa-file"></i> {}</span>', value, filename)

                    url, label, full_path = self._get_gcp_console_link(value, is_file=not is_dir)
                    if not url:
                        return value

                    icon = "fa-folder" if is_dir else "fa-file"
                    return format_html(
                        '<a href="{}" target="_blank" title="{}"><i class="fas {}"></i> {}</a>',
                        url,
                        full_path,
                        icon,
                        label,
                    )

                display_func.short_description = field_name.replace("_", " ").title()
                display_func.admin_order_field = field_name
                return display_func

            setattr(cls, method_name, make_display_func())

            if method_name not in cls.list_display:
                cls.list_display.append(method_name)

            cls._gcs_view_fields_meta[name] = is_directory

    def _get_gcp_console_link(
        self, path: str, bucket_name: str | None = None, is_file: bool = False
    ) -> tuple[str, str, str]:
        """
        Generate a GCP Console link for a GCS path.

        Create a link to the Google Cloud Console for viewing a file or directory in a GCS bucket.

        :param path: The path to the file or directory.
        :param bucket_name: The name of the bucket, if not included in the path.
        :param is_file: Whether the path points to a file (True) or directory (False).

        :raise: ValueError: If the path is in gs:// format but malformed.

        :return: A tuple of (url, label, full_path).
        """
        if not path:
            return "", "", path

        if path.startswith("gs://"):
            try:
                bucket, file_path = gcp.get_bucket_name_and_file_path_from_gc_path(path)
                full_path = path
            except ValueError:
                return "", "", path
        else:
            bucket = bucket_name or settings.BUCKET_NAME_PRIVATE
            file_path = path.lstrip("/")
            full_path = f"gs://{bucket}/{file_path}"

        if is_file:
            url = f"https://console.cloud.google.com/storage/browser/_details/{bucket}/{file_path}"
        else:
            url = f"https://console.cloud.google.com/storage/browser/{bucket}/{file_path}"

        return url, "View in Bucket", full_path

    def formfield_for_dbfield(self, db_field: Any, request: Any, **kwargs) -> Any:
        """
        Customize form fields for GCS viewable fields.

        Override the default formfield_for_dbfield method to add GCS view buttons to form fields
        that are in the list_gcs_viewable list.

        :param db_field: The database field being rendered.
        :param request: The current request object.
        :param kwargs: Additional keyword arguments for the form field.

        :return: The customized form field.
        """
        formfield = super().formfield_for_dbfield(db_field, request, **kwargs)

        # Inject GCS view button into field help text if it's in list_gcs_viewable
        if hasattr(self, "_gcs_view_fields_meta") and db_field.name in self._gcs_view_fields_meta:
            is_dir = self._gcs_view_fields_meta[db_field.name]
            object_id = request.resolver_match.kwargs.get("object_id")
            instance = self.model.objects.filter(pk=object_id).first() if object_id else None
            value = getattr(instance, db_field.name, None) if instance else None

            if value:
                url, label, full_path = self._get_gcp_console_link(str(value), is_file=not is_dir)
                if url:
                    html = f'<br><a href="{url}" target="_blank" class="button" style="margin-top:4px;">🔗 {label}</a>'
                    formfield.help_text = mark_safe((formfield.help_text or "") + html)

        return formfield
