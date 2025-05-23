"""
Admin module for cell information management.
"""

import logging
from typing import Any
from urllib.parse import parse_qs, urlparse

from django.contrib import admin, messages
from django.http import HttpRequest, HttpResponse
from django.shortcuts import redirect, render
from django.urls import reverse
from django.utils.http import urlencode
from django.utils.safestring import mark_safe
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin
from unfold.contrib.filters.admin import RangeNumericFilter, RelatedDropdownFilter
from unfold.decorators import action

from cellarium.nexus.backend.cell_management.admin import constants, filters, forms
from cellarium.nexus.backend.cell_management.admin import utils as admin_utils
from cellarium.nexus.backend.cell_management.models import BigQueryDataset, CellInfo
from cellarium.nexus.backend.core.admin.helpers.change_lists import BigQueryCountPaginator

logger = logging.getLogger(__name__)


@admin.register(CellInfo)
class CellInfoAdmin(ModelAdmin):
    """
    Admin interface for managing cell information.

    Provides functionality to filter and search cell data, including numeric range filters
    for fields like id and total_mrna_umis. Prevents direct creation and editing of cell info
    instances.
    """

    list_display = (
        "id",
        "original_id",
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
    )
    list_filter_submit = True
    list_filter = (
        # Text filters for fields without ontology term IDs
        ("id", RangeNumericFilter),
        ("total_mrna_umis", RangeNumericFilter),
        filters.SexDropdownFilter,
        filters.DonorDropdownFilter,
        filters.TagDropdownFilter,
        filters.SuspensionTypeDropdownFilter,
        filters.OrganismDropdownFilter,
        filters.CellTypeDropdownFilter,
        filters.DiseaseDropdownFilter,
        # Ontology term ID filters with specialized filter
        ("cell_type_ontology_term_id", filters.OntologyTermFilter),
        ("assay_ontology_term_id", filters.OntologyTermFilter),
        ("development_stage_ontology_term_id", filters.OntologyTermFilter),
        ("tissue_ontology_term_id", filters.OntologyTermFilter),
        ("disease_ontology_term_id", filters.OntologyTermFilter),
        ("organism_ontology_term_id", filters.OntologyTermFilter),
        ("self_reported_ethnicity_ontology_term_id", filters.OntologyTermFilter),
        # Related filters
        ("bigquery_dataset", RelatedDropdownFilter),
    )
    ordering = ("-id",)
    readonly_fields = ("id",)
    actions_list = ["extract_curriculum_action"]

    def has_add_permission(self, request: HttpRequest) -> bool:
        """
        Disable the ability to create new cell info instances directly.

        :param request: The HTTP request

        :return: False to prevent direct creation
        """
        return False

    def has_change_permission(self, request: HttpRequest, obj=None) -> bool:
        """
        Disable the ability to edit cell info instances.

        :param request: The HTTP request
        :param obj: The object being changed

        :return: False to prevent editing
        """
        return False

    def changelist_view(self, request, extra_context=None):
        if "bigquery_dataset__id__exact" not in request.GET:
            bigquery_dataset_first = BigQueryDataset.objects.first()
            base_url = reverse(f"admin:{self.model._meta.app_label}_{self.model._meta.model_name}_changelist")
            query_string = urlencode({"bigquery_dataset__id__exact": bigquery_dataset_first.id})
            return redirect(f"{base_url}?{query_string}")
        return super().changelist_view(request, extra_context=extra_context)

    def get_queryset(self, request: HttpRequest):
        return self.model.objects.none()

    @staticmethod
    def _get_filters_and_bigquery_dataset(request: HttpRequest) -> tuple[dict[str, Any], BigQueryDataset | None]:
        # Start with any filters in current GET
        query_params: dict[str, list[str]] = {k: request.GET.getlist(k) for k in request.GET}

        # If empty (e.g. follow-up view from action), try referrer
        if not query_params:
            referer = request.META.get("HTTP_REFERER", "")
            if referer:
                query_params = parse_qs(urlparse(referer).query)

        # Update request.GET to include filters so downstream logic works
        request.GET = request.GET.copy()
        request.GET.update(query_params)

        bq_filters, bigquery_dataset = admin_utils.extract_filters_from_django_admin_request(request=request)
        return bq_filters, bigquery_dataset

    def get_paginator(self, request, queryset, per_page, orphans=0, allow_empty_first_page=True):
        bq_filters, bigquery_dataset = self._get_filters_and_bigquery_dataset(request)
        return BigQueryCountPaginator(
            object_list=queryset,
            per_page=per_page,
            orphans=orphans,
            allow_empty_first_page=allow_empty_first_page,
            bq_filters=bq_filters,
            bigquery_dataset_name=bigquery_dataset.name if bigquery_dataset else None,
        )

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
        query_params = parse_qs(urlparse(referer).query)

        original_filters = {k: v for k, v in query_params.items()}

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
            extract_bin_keys = form.cleaned_data["extract_bin_keys"]

            pipeline_url = admin_utils.submit_extract_pipeline(
                feature_schema=feature_schema,
                name=name,
                extract_bin_size=extract_bin_size,
                categorical_column_count_limit=categorical_column_count_limit,
                extract_bin_keys=extract_bin_keys,
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
