"""
Admin module for cell information management.
"""

import logging
from urllib.parse import parse_qs, urlparse

from django.contrib import admin, messages
from django.http import HttpRequest, HttpResponse
from django.shortcuts import redirect, render
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter
from unfold.decorators import action

from cellarium.nexus.backend.cell_management.admin import constants, forms
from cellarium.nexus.backend.cell_management.admin import utils as admin_utils
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
    actions_list = ["extract_curriculum_action"]

    @action(description=_("Extract Curriculum"), url_path="extract-curriculum")
    def extract_curriculum_action(self, request: HttpRequest) -> HttpResponse:
        """
        Submit extract curriculum job

        :param request: The HTTP request

        :raise ValidationError: If form validation fails
        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

        :return: HTTP response
        """
        referer = request.META.get("HTTP_REFERER", "")
        # if referer and "?" in referer:
        query_params = parse_qs(urlparse(referer).query)
        original_filters = {k: v[0] for k, v in query_params.items()}
        logger.info(f"Original filter parameters: {original_filters}")

        request.GET = request.GET.copy()
        request.GET.update(original_filters)

        filters, bigquery_dataset = admin_utils.extract_filters_from_django_admin_request(request=request)

        # Create initial form data with extracted filters and dataset
        initial_data = {"filters": filters, "bigquery_dataset": bigquery_dataset}

        form = forms.ExtractCurriculumForm(request.POST or None, initial=initial_data)

        if request.method == "POST" and form.is_valid():
            # Get form data
            feature_schema = form.cleaned_data["feature_schema"]
            name = form.cleaned_data["name"]
            extract_bin_size = form.cleaned_data["extract_bin_size"]
            metadata_extra_columns = form.cleaned_data["metadata_extra_columns"]
            filters = form.cleaned_data["filters"] or {}
            bigquery_dataset = form.cleaned_data["bigquery_dataset"]
            categorical_column_count_limit = form.cleaned_data["categorical_column_count_limit"]

            pipeline_url = admin_utils.submit_extract_pipeline(
                feature_schema=feature_schema,
                name=name,
                extract_bin_size=extract_bin_size,
                categorical_column_count_limit=categorical_column_count_limit,
                filters=filters,
                bigquery_dataset=bigquery_dataset,
                creator_id=request.user.id,
                metadata_extra_columns=metadata_extra_columns,
            )

            messages.success(
                request=request, message=mark_safe(constants.EXTRACT_PIPELINE_SUCCESS_MESSAGE.format(pipeline_url))
            )
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
