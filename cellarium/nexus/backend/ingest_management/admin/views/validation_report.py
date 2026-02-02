from pathlib import Path

import pandas as pd
from django.conf import settings
from django.contrib import admin, messages
from django.core.exceptions import ValidationError
from django.http import HttpRequest, HttpResponse
from django.shortcuts import redirect, render
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin, TabularInline
from unfold.decorators import action

from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.admin import constants, forms
from cellarium.nexus.backend.ingest_management.utils.export_csv import export_model_queryset_to_csv
from cellarium.nexus.backend.ingest_management.utils.workflow_utils_local import run_validate_and_sanitize


class ValidationReportItemInline(TabularInline):
    """
    Inline admin for ValidationReportItem.

    Allows viewing validation report items directly within
    the ValidationReport admin interface. All fields are read-only.
    """

    model = models.ValidationReportItem
    extra = 0
    readonly_fields = ("truncated_gcs_path", "validator_name", "is_valid", "message", "created_at")
    fields = ("truncated_gcs_path", "validator_name", "is_valid", "message", "created_at")
    can_delete = False
    max_num = 0
    show_change_link = True

    def truncated_gcs_path(self, obj):
        """
        Display a truncated version of the GCS path with a tooltip showing the full path.

        :param obj: ValidationReportItem instance

        :return: Truncated GCS path with tooltip
        """
        if not obj.input_file_gcs_path:
            return "-"

        # Create a truncated version with ellipsis in the middle
        max_length = 50
        if len(obj.input_file_gcs_path) > max_length:
            prefix = obj.input_file_gcs_path[:20]
            suffix = obj.input_file_gcs_path[-25:]
            truncated = f"{prefix}...{suffix}"
        else:
            truncated = obj.input_file_gcs_path

        # Return HTML with tooltip showing full path
        return mark_safe(f'<span title="{obj.input_file_gcs_path}">{truncated}</span>')

    truncated_gcs_path.short_description = _("Input File GCS Path")


@admin.register(models.ValidationReport)
class ValidationReportAdmin(ModelAdmin):
    """
    Admin interface for managing validation reports.

    Provides functionality to view and manage validation reports with
    their associated validation items.
    """

    list_display = ("id", "created_at", "creator", "item_count", "download_report_button")
    list_filter = ("created_at", "creator")
    search_fields = ("id", "creator__username")
    readonly_fields = ("created_at",)
    inlines = [ValidationReportItemInline]
    actions_list = ["validate_new_data"]
    actions_detail = ["export_items_as_csv"]

    def item_count(self, obj):
        """
        Return the count of validation items in this report.

        :param obj: The ValidationReport instance

        :return: Number of validation items
        """
        return obj.items.count()

    item_count.short_description = "Number of Items"

    def download_report_button(self, obj: models.ValidationReport):
        """
        Generate a download button for the validation report items.

        :param obj: The ValidationReport instance

        :return: HTML for a download button
        """
        url = f"/admin/ingest_management/validationreport/{obj.id}/export-items-as-csv/?download=true"
        return mark_safe(f'<a class="button" href="{url}">Download CSV</a>')

    download_report_button.short_description = ""  # No column header

    @action(description=_("Export Items as CSV"), url_path="export-items-as-csv")
    def export_items_as_csv(self, request: HttpRequest, object_id: str) -> HttpResponse:
        """
        Export all validation report items as a CSV file.

        Generate a CSV file containing all validation report items for the specified report.
        The CSV file is named using the validation report ID and current timestamp.

        :param request: The HTTP request
        :param object_id: The validation report ID

        :raise ValueError: If validation report is not found or if queryset is empty
        :raise AttributeError: If a field cannot be accessed on a model instance

        :return: HTTP response with CSV file attachment or redirect to the detail page
        """
        try:
            validation_report = models.ValidationReport.objects.get(id=object_id)
        except models.ValidationReport.DoesNotExist:
            messages.error(request, _("Validation report not found"))
            return redirect("admin:ingest_management_validationreport_changelist")

        # Check if this is a download request or just a page refresh after download
        if request.GET.get("download") == "true":
            # Get all validation report items for this report
            items = models.ValidationReportItem.objects.filter(report=validation_report)

            # Use the utility function to generate the CSV response
            # Exclude the report field as it's not needed in the export
            response = export_model_queryset_to_csv(
                queryset=items,
                filename_prefix="validation_report",
                identifier=str(object_id),
                exclude_fields=["report"],
            )
            return response
        else:
            # Set a success message and redirect back to the detail page
            messages.success(request, _("CSV export has been downloaded"))
            return redirect("admin:ingest_management_validationreport_change", object_id=object_id)

    @action(description=_("Validate New Data"), url_path="validate-new-data")
    def validate_new_data(self, request: HttpRequest) -> HttpResponse:
        """
        Validate and sanitize SOMA data.

        :param request: The HTTP request

        :raises ValidationError: If validation fails
        :raises ValueError: If dataset or schema not found
        :raises IOError: If there's an error with GCS operations

        :return: HTTP response
        """
        form = forms.ValidateNewDataChangeListActionForm(request.POST or None, request.FILES or None)

        if request.method == "POST" and form.is_valid():
            dataset_name = form.cleaned_data["omics_dataset"].name
            csv_file = form.cleaned_data["input_csv_file"]
            output_directory = form.cleaned_data["output_directory_uri"].rstrip("/")

            # Read CSV file to extract input paths
            df = pd.read_csv(csv_file)

            # Get first column (assuming GCS paths are in first column)
            input_paths = df.iloc[:, 0].tolist()

            if not input_paths:
                raise ValidationError(_("CSV file contains no data"))

            # Generate output paths from directory and input filenames
            output_paths = [f"{output_directory}/{Path(path).name}" for path in input_paths]

            # Create validation report
            validation_report = models.ValidationReport.objects.create(creator=request.user)

            # Call validation pipeline
            run_validate_and_sanitize(
                dataset_name=dataset_name,
                input_h5ad_uris=input_paths,
                output_h5ad_uris=output_paths,
                validation_report_id=validation_report.id,
                nexus_backend_api_url=settings.SITE_URL,
            )

            messages.success(
                request=request,
                message=_(
                    "SOMA validation pipeline completed successfully. " f"Validation Report ID: {validation_report.id}"
                ),
            )
            return redirect("admin:ingest_management_validationreport_change", object_id=validation_report.id)

        return render(
            request=request,
            template_name=constants.CHANGELIST_ACTION_FORM,
            context={
                "form": form,
                "title": _("Validate New Data"),
                "submit_button_title": _("Validate"),
                **self.admin_site.each_context(request),
            },
        )


@admin.register(models.ValidationReportItem)
class ValidationReportItemAdmin(ModelAdmin):
    """
    Admin interface for managing individual validation report items.

    Provides detailed view of validation results for specific files.
    """

    list_display = ("id", "report", "validator_name", "is_valid", "truncated_gcs_path", "created_at")
    list_filter = ("is_valid", "validator_name", "created_at", "report")
    search_fields = ("input_file_gcs_path", "validator_name", "message")
    readonly_fields = ("created_at", "input_file_gcs_path", "truncated_gcs_path")

    def truncated_gcs_path(self, obj):
        """
        Display a truncated version of the GCS path with a tooltip showing the full path.

        :param obj: ValidationReportItem instance

        :return: Truncated GCS path with tooltip
        """
        if not obj.input_file_gcs_path:
            return "-"

        # Create a truncated version with ellipsis in the middle
        max_length = 50
        if len(obj.input_file_gcs_path) > max_length:
            prefix = obj.input_file_gcs_path[:20]
            suffix = obj.input_file_gcs_path[-25:]
            truncated = f"{prefix}...{suffix}"
        else:
            truncated = obj.input_file_gcs_path

        # Return HTML with tooltip showing full path
        return mark_safe(f'<span title="{obj.input_file_gcs_path}">{truncated}</span>')

    truncated_gcs_path.short_description = _("Input File GCS Path")
