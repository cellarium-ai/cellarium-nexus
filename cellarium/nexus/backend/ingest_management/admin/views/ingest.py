import pandas as pd
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
from cellarium.nexus.backend.ingest_management.utils import column_mapping_utils, workflows_utils


class ObsColumnMappingInline(TabularInline):
    model = models.ObsColumnMapping
    extra = 0
    verbose_name = _("Obs Column Mapping")
    verbose_name_plural = _("Obs Column Mappings")
    fields = ("input_column", "schema_column")


class VarColumnMappingInline(TabularInline):
    model = models.VarColumnMapping
    extra = 0
    verbose_name = _("Var Column Mapping")
    verbose_name_plural = _("Var Column Mappings")
    fields = ("input_column", "schema_column")


@admin.register(models.ColumnMapping)
class ColumnMappingAdmin(ModelAdmin):
    list_display = ("name", "description", "created_at", "updated_at")
    search_fields = ("name", "description")
    readonly_fields = ("created_at", "updated_at")
    inlines = [ObsColumnMappingInline, VarColumnMappingInline]
    fieldsets = (
        (None, {"fields": ("name", "description")}),
        (_("Timestamps"), {"fields": ("created_at", "updated_at"), "classes": ("collapse",)}),
    )


@admin.register(models.IngestInfo)
class IngestInfoAdmin(ModelAdmin):
    """
    Admin interface for managing ingest file information.

    Provides functionality to track and manage file ingestion metadata and data ingestion actions.
    Prevents direct creation and editing of ingests - they can only be created through the
    'Ingest New Data' action.
    """

    # Prevent adding new instances through the admin interface
    def has_add_permission(self, request: HttpRequest) -> bool:
        """
        Disable the ability to create new ingest instances directly.

        :param request: The HTTP request

        :return: False to prevent direct creation
        """
        return False

    def has_change_permission(self, request: HttpRequest, obj=None) -> bool:
        """
        Disable the ability to edit ingest instances.

        :param request: The HTTP request
        :param obj: The object being changed

        :return: False to prevent editing
        """
        return False

    list_display = (
        "id",
        "nexus_uuid",
        "status",
        "omics_dataset",
        "ingest_start_timestamp",
        "ingest_finish_timestamp",
    )
    search_fields = ("nexus_uuid", "status")
    list_filter = ("status", "omics_dataset")
    ordering = ("-ingest_start_timestamp",)
    readonly_fields = ("ingest_start_timestamp", "ingest_finish_timestamp", "nexus_uuid")
    fieldsets = (
        (
            None,
            {
                "fields": (
                    "nexus_uuid",
                    "status",
                    "omics_dataset",
                    "metadata_extra",
                    "ingest_start_timestamp",
                    "ingest_finish_timestamp",
                )
            },
        ),
    )
    actions_list = ["ingest_new_data", "validate_new_data"]

    @staticmethod
    def __get_validation_methods(gencode_version: str | None) -> list[str]:
        validation_methods = ["nexus.omics_datastore.bq_ops.validate.validate_raw_counts"]

        if gencode_version is not None:
            gencode_validation_method = (
                f"nexus.omics_datastore.bq_ops.validate.cellarium_validate_gencode_{gencode_version}"
            )
            validation_methods.append(gencode_validation_method)

        return validation_methods

    @action(description=_("Ingest New Data"), url_path="ingest-new-data")
    def ingest_new_data(self, request: HttpRequest) -> HttpResponse:
        """
        Ingest new data into the system using Kubeflow pipeline.

        :param request: The HTTP request

        :raise ValidationError: If ingestion fails
        :raise IOError: If there's an error writing configs to GCS

        :return: HTTP response
        """
        form = forms.IngestNewDataChangeListActionForm(request.POST or None, request.FILES or None)

        if request.method == "POST" and form.is_valid():
            csv_file = form.cleaned_data["ingest_csv_file"]
            omics_dataset = form.cleaned_data["omics_dataset"]
            column_mapping_obj = form.cleaned_data.get("column_mapping")
            # Get the selected genecode version
            gencode_version = form.cleaned_data["gencode_version"]

            validation_methods = self.__get_validation_methods(gencode_version=gencode_version)
            df = pd.read_csv(csv_file)
            if not all(col in df.columns for col in constants.REQUIRED_CSV_FILE_COLUMNS):
                raise ValidationError(f"CSV must contain columns: {', '.join(constants.REQUIRED_CSV_FILE_COLUMNS)}")

            column_mapping = column_mapping_utils.create_column_mapping(column_mapping_obj=column_mapping_obj)
            pipeline_url = workflows_utils.submit_ingest_pipeline(
                df_ingest_file_info=df,
                column_mapping=column_mapping,
                omics_dataset=omics_dataset,
                validation_methods=validation_methods,
            )

            messages.success(
                request=request, message=mark_safe(constants.INGEST_PIPELINE_SUCCESS_MESSAGE.format(pipeline_url))
            )
            return redirect("admin:ingest_management_ingestinfo_changelist")

        return render(
            request=request,
            template_name=constants.CHANGELIST_ACTION_FORM,
            context={
                "form": form,
                "title": _("Ingest New Data"),
                "submit_button_title": _("Ingest"),
                **self.admin_site.each_context(request),
            },
        )

    @action(description=_("Validate New Data"), url_path="validate-new-data")
    def validate_new_data(self, request: HttpRequest) -> HttpResponse:
        """
        Redirect to the ValidationReportAdmin's validate_new_data action.

        :param request: The HTTP request

        :return: HTTP response with redirect to ValidationReportAdmin action
        """
        return redirect("admin:ingest_management_validationreport_validate_new_data")
