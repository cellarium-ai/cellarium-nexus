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
from cellarium.nexus.backend.ingest_management.admin import utils as admin_utils

# Constants at the top


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
    """

    list_display = (
        "id",
        "nexus_uuid",
        "status",
        "bigquery_dataset",
        "ingest_start_timestamp",
        "ingest_finish_timestamp",
    )
    search_fields = ("nexus_uuid", "status")
    list_filter = ("status", "bigquery_dataset")
    ordering = ("-ingest_start_timestamp",)
    readonly_fields = ("ingest_start_timestamp", "ingest_finish_timestamp", "nexus_uuid")
    fieldsets = (
        (
            None,
            {
                "fields": (
                    "nexus_uuid",
                    "status",
                    "bigquery_dataset",
                    "metadata_extra",
                    "ingest_start_timestamp",
                    "ingest_finish_timestamp",
                )
            },
        ),
    )
    actions_list = ["ingest_new_data"]

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
            bigquery_dataset = form.cleaned_data["bigquery_dataset"]
            column_mapping_obj = form.cleaned_data.get("column_mapping")

            df = pd.read_csv(csv_file)
            if not all(col in df.columns for col in constants.REQUIRED_CSV_FILE_COLUMNS):
                raise ValidationError(f"CSV must contain columns: {', '.join(constants.REQUIRED_CSV_FILE_COLUMNS)}")

            column_mapping = admin_utils.create_column_mapping(column_mapping_obj=column_mapping_obj)
            pipeline_url = admin_utils.submit_ingest_pipeline(
                df_ingest_file_info=df, column_mapping=column_mapping, bigquery_dataset=bigquery_dataset
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
