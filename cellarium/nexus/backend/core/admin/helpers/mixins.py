"""
Admin mixins for cell management.
"""

from collections import Counter
from typing import Any

from django.contrib.admin import helpers
from django.contrib.admin.utils import model_ngettext
from django.core.exceptions import PermissionDenied
from django.db import models
from django.http import HttpRequest, HttpResponse
from django.template.response import TemplateResponse
from django.utils.translation import gettext


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
                    self.log_deletion(request=request, obj=obj, obj_display=obj_display)
                queryset.delete()
                self.message_user(
                    request=request,
                    message=gettext("Successfully deleted %(count)d %(items)s.")
                    % {"count": n, "items": model_ngettext(opts=self.opts, n=n)},
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
