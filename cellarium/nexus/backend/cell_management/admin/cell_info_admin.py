"""
Admin module for cell information management.
"""

import logging
from urllib.parse import parse_qs, urlparse

from django.contrib import admin, messages
from django.http import HttpRequest, HttpResponse
from django.shortcuts import redirect, render
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter
from unfold.decorators import action

from cellarium.nexus.backend.cell_management.admin import constants, forms, utils as admin_utils

from cellarium.nexus.backend.cell_management.admin.filters import MultiValueTextFilter, OntologyTermFilter
from cellarium.nexus.backend.cell_management.models import BigQueryDataset, CellInfo

logger = logging.getLogger(__name__)


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
    actions_list = ["extract_data_action"]

    @action(description=_("Extract Data"), url_path="extract-data")
    def extract_data_action(self, request: HttpRequest) -> HttpResponse:
        """
        Prepare extract tables for data extraction.

        :param request: The HTTP request

        :raise ValidationError: If form validation fails
        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

        :return: HTTP response
        """
        filters, bigquery_dataset = admin_utils.extract_filters_from_django_admin_request(request=request)

        # If no dataset from filters, try to get the default one
        if not bigquery_dataset:
            bigquery_dataset = BigQueryDataset.objects.get_default_dataset()
            if not bigquery_dataset:
                dataset_count = BigQueryDataset.objects.count()
                if dataset_count == 0:
                    messages.error(request, constants.NO_DATASETS_ERROR)
                    return redirect("admin:cell_management_cellinfo_changelist")
                else:
                    messages.error(request, _(constants.MULTIPLE_DATASETS_ERROR))
                    return redirect("admin:cell_management_cellinfo_changelist")

        referer = request.META.get("HTTP_REFERER", "")
        if referer and "?" in referer:
            query_params = parse_qs(urlparse(referer).query)
            original_filters = {k: v[0] for k, v in query_params.items()}
            logger.info(f"Original filter parameters: {original_filters}")

            request.GET = request.GET.copy()
            request.GET.update(original_filters)

            filters, bq_dataset_unused = admin_utils.extract_filters_from_django_admin_request(request=request)

        # Create initial form data with extracted filters and dataset
        initial_data = {"filters": filters, "bigquery_dataset": bigquery_dataset}

        form = forms.PrepareExtractTablesForm(request.POST or None, initial=initial_data)

        if request.method == "POST" and form.is_valid():
            # Get form data
            feature_schema = form.cleaned_data["feature_schema"]
            extract_table_prefix = form.cleaned_data["extract_table_prefix"]
            extract_bin_size = form.cleaned_data["extract_bin_size"]
            filters = form.cleaned_data["filters"] or {}

            admin_utils.submit_extract_pipeline(
                feature_schema=feature_schema,
                extract_table_prefix=extract_table_prefix,
                extract_bin_size=extract_bin_size,
                filters=filters,
                bigquery_dataset=bigquery_dataset,
            )

            messages.success(request=request, message=constants.EXTRACT_SUCCESS_MESSAGE)
            return redirect("admin:cell_management_cellinfo_changelist")

        return render(
            request=request,
            template_name=constants.CHANGELIST_ACTION_FORM,
            context={
                "form": form,
                "title": constants.PREPARE_EXTRACT_TABLES_TITLE,
                "submit_button_title": constants.PREPARE_BUTTON_TITLE,
                **self.admin_site.each_context(request),
            },
        )
