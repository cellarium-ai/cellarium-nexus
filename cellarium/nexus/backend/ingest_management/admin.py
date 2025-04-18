import datetime
import os
import secrets

from django.conf import settings
from django.contrib import admin, messages
from django.core.exceptions import ValidationError
from django.http import HttpRequest, HttpResponse
from django.shortcuts import redirect, render
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin, TabularInline
from unfold.decorators import action
import pandas as pd

from cellarium.nexus.shared import utils
from cellarium.nexus.workflows.kubeflow.component_configs import IngestTaskConfig
from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline
from cellarium.nexus.workflows.kubeflow.pipelines import ingest_data_pipeline
from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.forms import IngestNewDataChangeListActionForm


# Constants at the top
BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT = (
    '<a href="{url_link}" target="_blank" style="text-decoration: underline;">View in BigQuery</a>'
)
BIGQUERY_SUCCESS_MESSAGE_TEXT = _("BigQuery dataset in GCP was created successfully.")
CHANGELIST_ACTION_FORM = "admin/custom_templates/changelist_action_with_form.html"
REQUIRED_CSV_FILE_COLUMNS = ["gcs_file_path"]

# Form titles and messages
PREPARE_EXTRACT_TABLES_TITLE = _("Prepare Extract Tables")
PREPARE_BUTTON_TITLE = _("Prepare")
NO_DATASETS_ERROR = _("No BigQuery datasets available")
MULTIPLE_DATASETS_ERROR = _("Multiple BigQuery datasets exist. Please select a BigQuery dataset for extraction.")
EXTRACT_SUCCESS_MESSAGE = _("Extract tables prepared successfully")


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
        form = IngestNewDataChangeListActionForm(request.POST or None, request.FILES or None)

        if request.method == "POST" and form.is_valid():
            csv_file = form.cleaned_data["ingest_csv_file"]
            bigquery_dataset = form.cleaned_data["bigquery_dataset"]
            column_mapping_obj = form.cleaned_data.get("column_mapping")

            column_mapping = self._create_column_mapping(column_mapping_obj)

            df = pd.read_csv(csv_file)
            if not all(col in df.columns for col in REQUIRED_CSV_FILE_COLUMNS):
                raise ValidationError(f"CSV must contain columns: {', '.join(REQUIRED_CSV_FILE_COLUMNS)}")

            timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
            base_stage_dir = f"{settings.BACKEND_PIPELINE_DIR}/data-ingests/{timestamp}_{secrets.token_hex(12)}"

            # Create list for combined task configs
            task_configs = []

            for i, row in df.iterrows():
                # Extract file name without extension from gcs_file_path
                file_path = row["gcs_file_path"]
                file_name = os.path.basename(file_path)
                file_name_without_ext = os.path.splitext(file_name)[0]

                # Create unique stage dir for this file
                stage_dir = f"{base_stage_dir}/{secrets.token_hex(8)}_{file_name_without_ext[:10]}"

                tag = row["tag"] if "tag" in df.columns else None

                # Create combined config for this task
                task_configs.append(
                    IngestTaskConfig(
                        project_id=settings.GCP_PROJECT_ID,
                        nexus_backend_api_url=settings.SITE_URL,
                        bigquery_dataset=bigquery_dataset.name,
                        data_source_path=file_path,
                        bucket_name=settings.BUCKET_NAME_PRIVATE,
                        ingest_bucket_path=stage_dir,
                        tag=tag,
                        metadata_columns=column_mapping,
                    )
                )

            # Save configs to GCS
            configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"

            task_config_paths = utils.workflows_configs.dump_configs_to_bucket(
                task_configs, configs_stage_dir
            )

            # Submit pipeline
            submit_pipeline(
                pipeline_component=ingest_data_pipeline,
                display_name=f"Nexus Ingest Data - {bigquery_dataset.name}",
                gcp_project=settings.GCP_PROJECT_ID,
                pipeline_kwargs={
                    "ingest_task_configs": task_config_paths,
                },
                service_account=settings.PIPELINE_SERVICE_ACCOUNT,
                pipeline_root_path=settings.PIPELINE_ROOT_PATH,
            )

            messages.success(request, _("Data ingestion pipeline started successfully"))
            return redirect("admin:ingest_management_ingestinfo_changelist")

        return render(
            request,
            CHANGELIST_ACTION_FORM,
            {
                "form": form,
                "title": _("Ingest New Data"),
                "submit_button_title": _("Ingest"),
                **self.admin_site.each_context(request),
            },
        )

    def _create_column_mapping(self, column_mapping_obj: models.ColumnMapping | None) -> dict | None:
        """
        Create a column mapping dictionary from a ColumnMapping model instance.

        :param column_mapping_obj: ColumnMapping model instance or None
        :return: Dictionary with obs and var mappings or None if no mappings exist
        """
        if not column_mapping_obj:
            return None

        result = {}

        obs_mapping = {m.input_column: m.schema_column for m in column_mapping_obj.obs_mappings.all()}
        if obs_mapping:
            result["obs_mapping"] = obs_mapping

        var_mapping = {m.input_column: m.schema_column for m in column_mapping_obj.var_mappings.all()}
        if var_mapping:
            result["var_mapping"] = var_mapping

        return result or None
