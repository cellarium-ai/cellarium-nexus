"""
Admin mixins for cell management.
"""

from collections import Counter
from typing import Any

from django.conf import settings
from django.contrib.admin import helpers
from django.contrib.admin.utils import model_ngettext
from django.core.exceptions import PermissionDenied
from django.db import models
from django.http import HttpRequest, HttpResponse
from django.template.response import TemplateResponse
from django.utils.html import format_html, mark_safe
from django.utils.translation import gettext

from cellarium.nexus.shared.utils import gcp


class CountRelatedObjectsDeleteMixin:
    """
    Mixin to customize the delete confirmation page to show counts of related objects
    instead of listing all items individually.

    This mixin overrides both the single object delete view and the bulk delete action
    to show a summary of related objects that will be deleted.
    """

    actions = ["delete_selected"]

    def delete_view(
        self, request: HttpRequest, object_id: str, extra_context: dict[str, Any] | None = None
    ) -> HttpResponse:
        """
        Override the default delete_view to show counts of related objects.

        :param request: The HTTP request
        :param object_id: The ID of the object to delete
        :param extra_context: Extra context to pass to the template

        :raise: Http404 if object does not exist

        :return: HTTP response
        """
        opts = self.model._meta
        obj = self.get_object(request=request, object_id=object_id)

        if obj is None:
            return self._get_obj_does_not_exist_redirect(request=request, opts=opts, object_id=object_id)

        # If the user has already confirmed the deletion, process it as normal
        if request.method == "POST":
            return super().delete_view(request=request, object_id=object_id, extra_context=extra_context)

        # Get related counts directly using SQL for better performance
        related_counts = self._get_related_counts(obj=obj)

        # Get the back URL from the referrer or fallback to the changelist
        back_url = request.META.get("HTTP_REFERER", "")
        if not back_url or "delete" in back_url:
            back_url = f"../"

        context = {
            **self.admin_site.each_context(request=request),
            "title": gettext("Are you sure?"),
            "object": obj,
            "object_name": str(opts.verbose_name),
            "opts": opts,
            "app_label": opts.app_label,
            "related_counts": dict(related_counts),
            "back_url": back_url,
            **(extra_context or {}),
        }

        return TemplateResponse(
            request=request,
            template="admin/custom_templates/count_related_objects_delete_confirmation.html",
            context=context,
        )

    def delete_selected(self, request: HttpRequest, queryset: models.QuerySet) -> HttpResponse | None:
        """
        Override the default delete_selected action to show counts of related objects.

        :param request: The HTTP request
        :param queryset: The queryset of objects to delete

        :raise: PermissionDenied if user does not have delete permission

        :return: HTTP response or None
        """
        opts = self.model._meta
        app_label = opts.app_label

        # Check that the user has delete permission
        if not self.has_delete_permission(request):
            raise PermissionDenied

        # If the user has already confirmed the deletion, process it as normal
        if request.POST.get("post") and request.method == "POST":
            n = queryset.count()
            if n:
                for obj in queryset:
                    obj_display = str(obj)
                    self.log_deletion(request=request, obj=obj, object_repr=obj_display)
                queryset.delete()
                self.message_user(
                    request=request,
                    message=gettext("Successfully deleted %(count)d %(items)s.")
                    % {"count": n, "items": model_ngettext(self.opts, n)},
                )
            return None

        # Collect all related objects for all selected objects
        all_related_counts = Counter()
        for obj in queryset:
            related_counts = self._get_related_counts(obj=obj)
            all_related_counts.update(related_counts)

        # The rest of this is based on ModelAdmin.delete_selected_confirmation
        objects_name = model_ngettext(obj=queryset)

        context = {
            **self.admin_site.each_context(request=request),
            "title": gettext("Are you sure?"),
            "queryset": queryset,
            "objects_name": objects_name,
            "opts": opts,
            "app_label": app_label,
            "action_checkbox_name": helpers.ACTION_CHECKBOX_NAME,
            "related_counts": dict(all_related_counts),
        }

        return TemplateResponse(
            request=request,
            template="admin/custom_templates/count_related_objects_delete_selected_confirmation.html",
            context=context,
        )

    # Set the short description for the action
    delete_selected.short_description = gettext("Delete selected %(verbose_name_plural)s")

    def _get_related_counts(self, obj: models.Model, processed_models: set | None = None) -> Counter:
        """
        Get counts of related objects that will be deleted when deleting the given object.
        This method recursively traverses related objects to get counts at each level.

        :param obj: The object to get related counts for
        :param processed_models: Set of models that have already been processed (to avoid infinite recursion)

        :raise: AttributeError if related object does not have expected attributes

        :return: Counter of model verbose names to counts
        """
        if processed_models is None:
            processed_models = set()

        model = obj.__class__
        model_meta = model._meta
        counts = Counter()

        # Skip if we have already processed this model to avoid infinite recursion
        if model in processed_models:
            return counts

        # Add this model to processed models
        processed_models = processed_models.union({model})

        # Get all related objects
        for related_object in model_meta.get_fields():
            # Skip if not a relation or a many-to-many relation
            if not hasattr(related_object, "related_model") or related_object.many_to_many:
                continue

            # Skip if the related model is the same as the current model (self-reference)
            if related_object.related_model == model:
                continue

            # Skip if we have already processed this model
            if related_object.related_model in processed_models:
                continue

            # Get the related manager
            if hasattr(related_object, "get_accessor_name"):
                # Reverse relation
                accessor_name = related_object.get_accessor_name()
                if not hasattr(obj, accessor_name):
                    continue
                related_manager = getattr(obj, accessor_name)
            else:
                # Forward relation (ForeignKey)
                continue

            # Count related objects
            count = related_manager.count()
            if count > 0:
                related_model_name = related_object.related_model._meta.verbose_name_plural
                counts[related_model_name] = count

                # For nested relationships, sample a few objects to estimate counts
                # This is more efficient than loading all objects
                sample_size = min(5, count)
                if sample_size > 0:
                    sample_objects = related_manager.all()[:sample_size]
                    nested_processed = processed_models.copy()

                    # Get counts for each sampled object and scale up
                    for sample_obj in sample_objects:
                        nested_counts = self._get_related_counts(obj=sample_obj, processed_models=nested_processed)
                        # Scale the counts based on the sample size
                        if sample_size < count:
                            scaling_factor = count / sample_size
                            for k, v in nested_counts.items():
                                nested_counts[k] = int(v * scaling_factor)
                        counts.update(nested_counts)

        return counts


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
