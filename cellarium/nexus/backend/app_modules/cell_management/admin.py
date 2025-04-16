import os
import csv
import logging
from typing import Sequence
import datetime

import pandas as pd
from django import forms
from django.conf import settings
from django.contrib import admin, messages
from django.core.exceptions import ValidationError
from django.db.models import QuerySet
from django.http import Http404, HttpRequest, HttpResponse, JsonResponse
from django.shortcuts import redirect, render
from django.urls import path, reverse
from django.utils.html import format_html, format_html_join
from django.utils.translation import gettext_lazy as _
from nexus.backend.app_modules.cell_management.forms import (
    CreateSchemaFromCSVForm,
    IngestNewDataChangeListActionForm,
    PrepareExtractTablesForm,
)
from nexus.backend.app_modules.cell_management.models import (
    BigQueryDataset,
    CellFeatureInfo,
    CellInfo,
    ColumnMapping,
    FeatureInfo,
    FeatureSchema,
    IngestInfo,
    ObsColumnMapping,
    VarColumnMapping,
)

from cellarium.nexus.workflows.kubeflow.component_configs import (
    CreateIngestFiles,
    IngestDataToBigQuery,
    BQOpsPrepareExtract,
    BQOpsExtract,
)

from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline
from cellarium.nexus.workflows.kubeflow.pipelines import ingest_data_pipeline, extract_data_pipeline

from nexus.backend.app_modules.cell_management.utils.custom_filters import MultiValueTextFilter, OntologyTermFilter
from nexus.backend.app_modules.cell_management.utils.filters import (
    extract_filters_from_django_admin_request,
    get_default_dataset,
)
from cellarium.nexus.shared import schemas, utils

# from nexus.omics_datastore.controller import NexusDataController
from nexus.omics_datastore.bq_ops import create_bq_tables
from unfold.admin import ModelAdmin, TabularInline
from unfold.contrib.filters.admin import (
    ChoicesDropdownFilter,
    RangeNumericFilter,
    RelatedDropdownFilter,
)
from unfold.decorators import action

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

logger = logging.getLogger(__name__)


@admin.register(BigQueryDataset)
class BigQueryDatasetAdmin(ModelAdmin):
    """
    Admin interface for managing BigQuery datasets.

    Provides functionality to create and manage BigQuery datasets in GCP.
    """

    list_display = ("id", "name", "description", "link_display")
    search_fields = ("name",)
    list_filter = ("name",)
    ordering = ("name",)
    readonly_fields = ("link",)
    fieldsets = ((None, {"fields": ("name", "description", "link")}),)

    def save_model(self, request: HttpRequest, obj: BigQueryDataset, form: forms.Form, change: bool) -> None:
        """
        Create a BigQuery dataset when a new record is created.

        :param request: The HTTP request
        :param obj: The BigQueryDataset instance being saved
        :param form: The form used to create/edit the instance
        :param change: Boolean indicating if this is a change to an existing record

        :raise Exception: If dataset creation fails
        """
        is_new = not change  # `change` is False if adding a new object

        if is_new:
            try:
                # Initialize the controller first
                # controller = NexusDataController(
                #     project_id=settings.GCP_PROJECT_ID,
                #     nexus_backend_api_url=settings.SITE_URL,
                #     bigquery_dataset=obj.name,
                # )

                # Create the BigQuery dataset and get the link
                # link_to_dataset = controller.create_bigquery_dataset(bigquery_dataset=obj.name, location="US")
                from google.cloud import bigquery

                link_to_dataset = create_bq_tables.create_bigquery_objects(
                    client=bigquery.Client(), project=settings.GCP_PROJECT_ID, dataset=obj.name, location="US"
                )
                # Assign the link to the object's 'link' field for display in admin
                obj.link = link_to_dataset
            except Exception as e:
                # Display an error message if dataset creation fails
                self.message_user(request, f"Failed to create BigQuery dataset: {e}", level=messages.ERROR)
                return  # Exit early to prevent saving the object without a valid link

        super().save_model(request=request, obj=obj, form=form, change=change)

        if is_new:
            html_link = format_html(BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT, url_link=obj.link)
            self.message_user(
                request=request,
                message=format_html(BIGQUERY_SUCCESS_MESSAGE_TEXT + " " + html_link),
                level=messages.SUCCESS,
            )

    def link_display(self, obj: BigQueryDataset) -> str:
        """
        Generate clickable link to BigQuery dataset.

        :param obj: The current BigQueryDataset object
        :return: HTML anchor tag with the link
        """
        if obj.link:
            return format_html(BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT, url_link=obj.link)
        return "-"

    link_display.short_description = "BigQuery Link"
    link_display.allow_tags = True  # Optional: For compatibility with older Django versions


@admin.register(CellInfo)
class CellInfoAdmin(ModelAdmin):
    """
    Admin interface for managing cell information.

    Provides functionality to filter and search cell data, including numeric range filters
    for fields like id and total_mrna_umis.
    """

    list_display = (
        "id",
        "original_id",
        "ingest",
        "donor_id",
        "cell_type",
        "assay",
        "development_stage",
        "tissue",
        "disease",
        "organism",
        "self_reported_ethnicity",
        "sex",
        "suspension_type",
        "total_mrna_umis",
        "tag",
        "ingest__bigquery_dataset",
    )
    list_filter_submit = True
    search_fields = ("id", "original_id", "cell_type", "assay", "organism", "tissue", "disease", "tag")
    list_filter = (
        # Text filters for fields without ontology term IDs
        ("total_mrna_umis", RangeNumericFilter),
        ("id", RangeNumericFilter),
        "sex",
        ("donor_id", MultiValueTextFilter),
        ("tag", MultiValueTextFilter),
        ("suspension_type", MultiValueTextFilter),
        # Ontology term ID filters with specialized filter
        ("cell_type_ontology_term_id", OntologyTermFilter),
        ("assay_ontology_term_id", OntologyTermFilter),
        ("development_stage_ontology_term_id", OntologyTermFilter),
        ("tissue_ontology_term_id", OntologyTermFilter),
        ("disease_ontology_term_id", OntologyTermFilter),
        ("organism_ontology_term_id", OntologyTermFilter),
        ("self_reported_ethnicity_ontology_term_id", OntologyTermFilter),
        # Related filters
        ("ingest__bigquery_dataset", RelatedDropdownFilter),
    )
    ordering = ("-id",)
    readonly_fields = ("id",)
    actions_list = ["extract_data_action", "debug_filters"]

    @action(description=_("Extract Data"), url_path="extract-data")
    def extract_data_action(self, request: HttpRequest) -> HttpResponse:
        """
        Prepare extract tables for data extraction.

        :param request: The HTTP request

        :raise ValidationError: If form validation fails
        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

        :return: HTTP response
        """
        # First try to get the BigQuery dataset
        bigquery_dataset = None

        # Try to get dataset from request parameters
        filters, bigquery_dataset = extract_filters_from_django_admin_request(request)

        # If no dataset from filters, try to get the default one
        if not bigquery_dataset:
            bigquery_dataset = get_default_dataset()
            if not bigquery_dataset:
                dataset_count = BigQueryDataset.objects.count()
                if dataset_count == 0:
                    messages.error(request, _(NO_DATASETS_ERROR))
                    return redirect("admin:cell_management_cellinfo_changelist")
                else:
                    messages.error(request, _(MULTIPLE_DATASETS_ERROR))
                    return redirect("admin:cell_management_cellinfo_changelist")

        # Get the original filter parameters from the referer URL if available
        original_filters = {}
        referer = request.META.get("HTTP_REFERER", "")
        if referer and "?" in referer:
            from urllib.parse import parse_qs, urlparse

            query_params = parse_qs(urlparse(referer).query)
            # Convert query params to a flat dictionary
            original_filters = {k: v[0] for k, v in query_params.items()}
            logger.info(f"Original filter parameters: {original_filters}")

            # Update request.GET with original filters
            # Create a mutable copy of request.GET
            request.GET = request.GET.copy()
            request.GET.update(original_filters)

            # Re-extract filters with updated request.GET
            filters, _ = extract_filters_from_django_admin_request(request)

        # Create initial form data with extracted filters and dataset
        initial_data = {
            "filters": filters,  # Pass the filters dict directly, not as JSON string
            "bigquery_dataset": bigquery_dataset,  # Pass the dataset object directly
        }

        form = PrepareExtractTablesForm(request.POST or None, initial=initial_data)

        if request.method == "POST" and form.is_valid():
            # Get form data
            feature_schema = form.cleaned_data["feature_schema"]
            extract_table_prefix = form.cleaned_data["extract_table_prefix"]
            filters = form.cleaned_data["filters"] or {}

            # Use the pre-selected dataset instead of getting it from the form
            # since we've already validated its existence

            # Convert FeatureSchema to schemas.FeatureSchema sequence
            features: Sequence[schemas.FeatureSchema] = [
                schemas.FeatureSchema(id=idx, symbol=feature.symbol, ensemble_id=feature.ensemble_id)
                for idx, feature in enumerate(feature_schema.features.all())
            ]

            # Hardcoded columns based on CellInfo schema
            obs_columns = [
                "donor_id",
                "cell_type",
                "assay",
                "development_stage",
                "tissue",
                "disease",
                "organism",
                "self_reported_ethnicity",
                "sex",
                "suspension_type",
                "total_mrna_umis",
                "cell_type_ontology_term_id",
                "assay_ontology_term_id",
                "development_stage_ontology_term_id",
                "tissue_ontology_term_id",
                "disease_ontology_term_id",
                "organism_ontology_term_id",
                "self_reported_ethnicity_ontology_term_id",
                "sex_ontology_term_id",
            ]

            # Construct extract bucket path
            extract_bucket_path = f"{settings.BACKEND_PIPELINE_DIR}/data_extracts/{extract_table_prefix}"

            # Create configs for the pipeline
            prepare_extract_config = BQOpsPrepareExtract(
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=bigquery_dataset.name,
                extract_table_prefix=extract_table_prefix,
                features=features,
                filters=filters,
                obs_columns=obs_columns,
                extract_bin_size=10000,
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                extract_bucket_path=extract_bucket_path,
            )

            extract_config = BQOpsExtract(
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=bigquery_dataset.name,
                extract_table_prefix=extract_table_prefix,
                bins=[0],  # Extract only bin 0 for now
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                extract_bucket_path=extract_bucket_path,
                obs_columns=obs_columns,
                max_workers=10,
            )

            # Save configs to GCS
            configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"

            prepare_extract_config_path = utils.workflows_configs.dump_configs_to_bucket(
                [prepare_extract_config], configs_stage_dir
            )[0]
            extract_config_paths = utils.workflows_configs.dump_configs_to_bucket([extract_config], configs_stage_dir)

            # Submit pipeline
            submit_pipeline(
                pipeline_component=extract_data_pipeline,
                display_name=f"Nexus Extract Data - {bigquery_dataset.name}",
                gcp_project=settings.GCP_PROJECT_ID,
                pipeline_kwargs={
                    "prepare_extract_config": prepare_extract_config_path,
                    "extract_configs": extract_config_paths,
                },
                service_account=settings.PIPELINE_SERVICE_ACCOUNT,
                pipeline_root_path=settings.PIPELINE_ROOT_PATH,
            )

            messages.success(request=request, message=EXTRACT_SUCCESS_MESSAGE)
            return redirect("admin:cell_management_cellinfo_changelist")

        return render(
            request,
            CHANGELIST_ACTION_FORM,
            {
                "form": form,
                "title": PREPARE_EXTRACT_TABLES_TITLE,
                "submit_button_title": PREPARE_BUTTON_TITLE,
                **self.admin_site.each_context(request),
            },
        )

    @action(description=_("Debug Filters"), url_path="debug-filters")
    def debug_filters(self, request: HttpRequest) -> HttpResponse:
        """
        Debug view to show filter parameters.

        :param request: The HTTP request
        :return: HTTP response with filter debug info
        """
        # Log all GET parameters
        logger.info(f"Debug Filters - GET params: {request.GET}")

        # Create a response with all GET parameters
        response_data = {"GET_params": dict(request.GET), "filter_params": {}}

        # Extract filter parameters for each ontology term filter
        for field_path in [
            "cell_type_ontology_term_id",
            "assay_ontology_term_id",
            "development_stage_ontology_term_id",
            "tissue_ontology_term_id",
            "disease_ontology_term_id",
            "organism_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
        ]:
            value = request.GET.get(field_path, "")
            exclude = request.GET.get(f"{field_path}_exclude", "") == "on"

            if value:
                # Check if it's a comma-separated list
                if "," in value:
                    values = [v.strip() for v in value.split(",") if v.strip()]
                    filter_kwargs = {f"{field_path}__in": values}
                else:
                    filter_kwargs = {f"{field_path}__exact": value}

                # Add to response data
                response_data["filter_params"][field_path] = {
                    "value": value,
                    "exclude": exclude,
                    "filter_kwargs": filter_kwargs,
                }

        # Return as JSON
        return JsonResponse(response_data)


@admin.register(IngestInfo)
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
            base_stage_dir = f"{settings.BACKEND_PIPELINE_DIR}/data-ingests"

            # Create ingest file configs for each row in the CSV
            create_ingest_configs = []
            stage_dirs = []
            for i, row in df.iterrows():
                # Extract file name without extension from gcs_file_path
                file_path = row["gcs_file_path"]
                file_name = os.path.basename(file_path)
                file_name_without_ext = os.path.splitext(file_name)[0]

                # Create unique stage dir for this file
                stage_dir = f"{base_stage_dir}/{timestamp}_{file_name_without_ext}"
                # stage_dirs.append(stage_dir)

                tag = row["tag"] if "tag" in df.columns else None
                create_ingest_configs.append(
                    CreateIngestFiles(
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

            # Create single ingest config for all files
            # Use base_stage_dir since IngestDataToBigQuery will look in all subdirs
            ingest_config = IngestDataToBigQuery(
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=bigquery_dataset.name,
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                ingest_bucket_path=base_stage_dir,
            )

            # Save configs to GCS
            configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"

            create_ingest_config_paths = utils.workflows_configs.dump_configs_to_bucket(
                create_ingest_configs, configs_stage_dir
            )
            # create_ingest_configs = [{"gcs_config_path": x} for x in create_ingest_config_paths]

            ingest_paths = utils.workflows_configs.dump_configs_to_bucket([ingest_config], configs_stage_dir)
            # Submit pipeline
            submit_pipeline(
                pipeline_component=ingest_data_pipeline,
                display_name=f"Nexus Ingest Data - {bigquery_dataset.name}",
                gcp_project=settings.GCP_PROJECT_ID,
                pipeline_kwargs={
                    "create_ingest_configs": create_ingest_config_paths,
                    "ingest_config": ingest_paths[0],
                },
                service_account=settings.PIPELINE_SERVICE_ACCOUNT,
                pipeline_root_path=settings.PIPELINE_ROOT_PATH,
            )

            messages.success(request, _("Data ingestion pipeline started successfully"))
            return redirect("admin:cell_management_ingestinfo_changelist")

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

    def _create_column_mapping(self, column_mapping_obj: ColumnMapping | None) -> dict | None:
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


@admin.register(CellFeatureInfo)
class CellFeatureInfoAdmin(ModelAdmin):
    """
    Admin interface for managing cell feature information.

    Provides functionality to filter and search cell feature data, including numeric range filters
    for the id field.
    """

    list_display = (
        "id",
        "ensemble_id",
        "symbol",
        "biotype",
        "is_filtered",
        "reference",
        "ingest",
        "tag",
    )
    search_fields = (
        "id",
        "ensemble_id",
        "symbol",
        "biotype",
        "reference",
        "tag",
    )
    list_filter = (
        ("id", RangeNumericFilter),
        ("biotype", MultiValueTextFilter),
        ("is_filtered", ChoicesDropdownFilter),
        ("tag", MultiValueTextFilter),
        ("ensemble_id", MultiValueTextFilter),
        ("symbol", MultiValueTextFilter),
        ("ingest", RelatedDropdownFilter),
    )
    ordering = ("-id",)
    readonly_fields = ("id",)


@admin.register(FeatureSchema)
class FeatureSchemaAdmin(ModelAdmin):
    """
    Admin interface for managing feature schemas.

    Provides functionality to create, download and manage feature schemas.
    """

    list_display = ("name", "feature_count", "get_download_button")
    search_fields = ("name",)
    ordering = ("name",)
    readonly_fields = ("features_display",)
    fields = ("name", "features_display")
    actions = ["download_features_csv_action"]
    actions_list = ("upload_features_csv",)

    def get_urls(self) -> list:
        """
        Add custom URLs for the admin interface.

        :return: List of URL patterns
        """
        urls = super().get_urls()
        info = self.model._meta.app_label, self.model._meta.model_name
        custom_urls = [
            path(
                "<path:object_id>/download-csv/",
                self.admin_site.admin_view(self.download_features_csv),
                name="%s_%s_download-csv" % info,
            ),
        ]
        return custom_urls + urls

    def feature_count(self, obj: FeatureSchema) -> int:
        """
        Get count of features in schema.

        :param obj: FeatureSchema instance
        :return: Count of features
        """
        return obj.features.count()

    feature_count.short_description = _("Feature Count")

    def get_download_button(self, obj: FeatureSchema) -> str:
        """
        Get download button for list view.

        :param obj: FeatureSchema instance
        :return: HTML for download button
        """
        info = self.model._meta.app_label, self.model._meta.model_name
        url = reverse(f"admin:{info[0]}_{info[1]}_download-csv", args=[obj.pk])
        return format_html('<a class="button" href="{}">{}</a>', url, _("Download CSV"))

    get_download_button.short_description = _("Download")

    def features_display(self, obj: FeatureSchema) -> str:
        """
        Display features with download button.

        :param obj: FeatureSchema instance
        :return: HTML for features display
        """
        features = obj.features.all().order_by("symbol")

        # Create a styled list of features
        feature_list = format_html(
            '<div style="margin-bottom: 15px;">'
            '<ul style="list-style-type: none; padding: 0; margin: 0;">{}</ul>'
            "</div>",
            format_html_join(
                "", '<li style="padding: 3px 0;">{} - {}</li>', ((f.ensemble_id, f.symbol) for f in features[:10])
            ),
        )

        if features.count() > 10:
            feature_list = format_html(
                '{}<div style="color: #666;">... and {} more features</div>', feature_list, features.count() - 10
            )

        # Add the download button with some spacing
        info = self.model._meta.app_label, self.model._meta.model_name
        download_url = reverse(f"admin:{info[0]}_{info[1]}_download-csv", args=[obj.pk])

        download_button = format_html(
            '<div style="margin-top: 10px;">'
            '<a href="{}" class="button" style="display: inline-block;">{}</a>'
            "</div>",
            download_url,
            _("Download CSV"),
        )

        return format_html("{}\n{}", feature_list, download_button)

    def download_features_csv(self, request: HttpRequest, object_id: int) -> HttpResponse:
        """
        Download features as CSV.

        :param request: The HTTP request
        :param object_id: Schema ID
        :return: CSV response
        """
        schema = self.get_object(request, object_id)
        if not schema:
            raise Http404(_("Schema does not exist"))

        features = schema.features.all().order_by("symbol")

        response = HttpResponse(content_type="text/csv")
        response["Content-Disposition"] = f'attachment; filename="{schema.name}_features.csv"'

        writer = csv.writer(response)
        writer.writerow(["ensemble_id", "symbol"])
        writer.writerows(features.values_list("ensemble_id", "symbol"))

        return response

    @admin.action(description=_("Download features as CSV"))
    def download_features_csv_action(self, request: HttpRequest, queryset: QuerySet) -> HttpResponse:
        """
        Admin action to download features as CSV.

        :param request: The HTTP request
        :param queryset: QuerySet of selected schemas
        :return: CSV file response
        """
        if len(queryset) != 1:
            self.message_user(request, _("Please select exactly one schema to download features from."), messages.ERROR)
            return

        schema = queryset.first()
        return self.download_features_csv(request, schema.pk)

    @action(description=_("Upload Features CSV"), url_path="upload-csv")
    def upload_features_csv(
        self, request: HttpRequest, permissions: tuple[str] = ("upload_features_csv",)
    ) -> HttpResponse:
        """
        Handle CSV upload for creating new schema.

        :param request: The HTTP request
        :param permissions: Required permissions

        :raise ValidationError: If CSV file is invalid

        :return: HTTP response
        """
        form = CreateSchemaFromCSVForm(data=request.POST or None, files=request.FILES or None)

        if request.method == "POST" and form.is_valid():
            try:
                schema_name = form.cleaned_data["schema_name"]
                csv_file = form.cleaned_data["csv_file"]

                df = pd.read_csv(csv_file)
                schema = FeatureSchema.objects.create(name=schema_name)

                # Get all existing features that match our CSV data
                existing_features = set(
                    FeatureInfo.objects.filter(
                        ensemble_id__in=df["ensemble_id"].tolist(), symbol__in=df["symbol"].tolist()
                    ).values_list("ensemble_id", "symbol")
                )

                # Create only features that don't exist yet
                new_features = []
                for index, row in df.iterrows():
                    if (row["ensemble_id"], row["symbol"]) not in existing_features:
                        new_features.append(FeatureInfo(ensemble_id=row["ensemble_id"], symbol=row["symbol"]))

                # Bulk create new features if any
                if new_features:
                    FeatureInfo.objects.bulk_create(new_features)

                # Get all features (both existing and newly created) for the schema
                all_features = FeatureInfo.objects.filter(
                    ensemble_id__in=df["ensemble_id"].tolist(), symbol__in=df["symbol"].tolist()
                )
                schema.features.add(*all_features)

                messages.success(
                    request,
                    _("Schema created successfully. Added {} features ({} new, {} existing)").format(
                        len(df), len(new_features), len(existing_features)
                    ),
                )
                return redirect("admin:cell_management_featureschema_changelist")
            except Exception as e:
                messages.error(request, _("Error creating schema: %s") % str(e))

        return render(
            request=request,
            template_name=CHANGELIST_ACTION_FORM,
            context={
                "form": form,
                "title": _("Upload Features CSV"),
                "submit_button_title": _("Upload"),
                **self.admin_site.each_context(request),
            },
        )


class ObsColumnMappingInline(TabularInline):
    model = ObsColumnMapping
    extra = 0
    verbose_name = _("Obs Column Mapping")
    verbose_name_plural = _("Obs Column Mappings")
    fields = ("input_column", "schema_column")


class VarColumnMappingInline(TabularInline):
    model = VarColumnMapping
    extra = 0
    verbose_name = _("Var Column Mapping")
    verbose_name_plural = _("Var Column Mappings")
    fields = ("input_column", "schema_column")


@admin.register(ColumnMapping)
class ColumnMappingAdmin(ModelAdmin):
    list_display = ("name", "description", "created_at", "updated_at")
    search_fields = ("name", "description")
    readonly_fields = ("created_at", "updated_at")
    inlines = [ObsColumnMappingInline, VarColumnMappingInline]
    fieldsets = (
        (None, {"fields": ("name", "description")}),
        (_("Timestamps"), {"fields": ("created_at", "updated_at"), "classes": ("collapse",)}),
    )

    def get_mapping_dict(self) -> dict:
        """Convert the relational mapping to the format expected by the ingest process."""
        obs_mapping = {m.input_column: m.schema_column for m in self.obs_mappings.all()}
        var_mapping = {m.input_column: m.schema_column for m in self.var_mappings.all()}
        return {"obs_mapping": obs_mapping, "var_mapping": var_mapping}
