"""
Admin module for cell information management.
"""

from django.contrib import admin
from django.contrib import messages
import django.forms as dj_forms
from django.http import JsonResponse as http_JsonResponse
import django.views.decorators.http as http_decorators
import django.views.generic as generic
import google.api_core as google_api_core
import json as py_json
import typing as t

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.admin.utils import bigquery_utils
from cellarium.nexus.backend.cell_management.admin import forms as admin_forms


class CellInfoAdminView(generic.TemplateView):
    """
    Render the Cell Info page within the admin site with BigQuery-backed counts.

    :param request: The HTTP request

    :raise: None

    :return: The rendered placeholder template using admin base
    """

    template_name = "admin/cell_management/cellinfo/cell_info_bigquery.html"

    def get_context_data(self, **kwargs):
        """
        Build template context with datasets from the database and BigQuery counts.

        :param kwargs: The keyword args passed by the TemplateView

        :raise: None

        :return: The template context dictionary
        """
        context = super().get_context_data(**kwargs)

        # Load datasets from the database
        bq_datasets_qs = cell_models.BigQueryDataset.objects.all().order_by("name")
        datasets: list[str] = [ds.name for ds in bq_datasets_qs]

        # Choose a default selection: use singleton default if present, else first
        default_ds_obj = bq_datasets_qs.get_default_dataset() if hasattr(bq_datasets_qs, "get_default_dataset") else None
        selected_dataset = default_ds_obj.name if default_ds_obj else (datasets[0] if datasets else "")

        # Compute counts with caching; handle BigQuery errors gracefully
        dataset_counts: dict[str, int] = {}
        manager = bigquery_utils.BigQueryCachedDataManager()
        try:
            for ds in datasets:
                dataset_counts[ds] = manager.get_cached_count_bq(dataset_name=ds, filters_dict={})
        except google_api_core.exceptions.GoogleAPICallError as exc:
            messages.error(
                request=self.request,
                message=(
                    "Unable to retrieve cell counts from BigQuery at the moment. "
                    "Please try again later or contact an administrator."
                ),
            )

        # Build server-rendered filters formset
        fields_meta = _get_cellinfo_filters_fields()
        field_choices: list[tuple[str, str]] = [
            (f["key"], f.get("label") or f["key"]) for f in fields_meta
        ]

        class _BaseFilterFormSet(dj_forms.formsets.BaseFormSet):
            def get_form_kwargs(self_inner, index: int) -> dict:
                return {
                    "field_choices": field_choices,
                }

        FilterRowFormSet = dj_forms.formset_factory(
            form=admin_forms.FilterRowForm,
            formset=_BaseFilterFormSet,
            extra=1,
            can_delete=True,
        )
        formset = FilterRowFormSet(prefix="filters")

        context.update({
            "datasets": datasets,
            "selected_dataset": selected_dataset,
            "dataset_counts": dataset_counts,
            "filters_formset": formset,
            "filters_fields_meta": fields_meta,
            **admin.site.each_context(self.request),
        })
        return context


def _get_cellinfo_filters_fields() -> list[dict]:
    """
    Return filters metadata used by the server-rendered formset and API endpoints.

    :raise: None

    :return: List of field metadata dictionaries
    """
    return [
        {
            "key": "cell_type",
            "label": "Cell Type",
            "type": "string",
            "operators": ["eq", "not_eq", "in", "not_in"],
            "suggest": True,
            "suggest_mode": "prefix",
            "suggest_min_chars": 2,
        },
        {
            "key": "assay",
            "label": "Assay",
            "type": "string",
            "operators": ["eq", "in"],
            "suggest": True,
            "suggest_mode": "prefix",
            "suggest_min_chars": 2,
        },
        {
            "key": "total_mrna_umis",
            "label": "Total mRNA UMIs",
            "type": "number",
            "operators": ["eq", "gt", "gte", "lt", "lte"],
            "suggest": False,
        },
        {
            "key": "is_high_quality",
            "label": "High Quality",
            "type": "boolean",
            "operators": ["eq", "not_eq"],
            "suggest": False,
        },
    ]


@http_decorators.require_GET
def cellinfo_filters_fields(request):
    """
    Return filters metadata for the filter builder UI.

    :param request: The HTTP request.

    :raise: None

    :return: A JSON response with fields metadata including operators and suggestion hints.
    """
    fields = _get_cellinfo_filters_fields()
    return http_JsonResponse(data={"fields": fields})


@http_decorators.require_POST
def cellinfo_filters_count(request):
    """
    Return a mocked count based on provided dataset and filters.

    :param request: The HTTP request.

    :raise ValueError: If payload is not valid JSON.

    :return: A JSON response containing a mocked count integer.
    """
    try:
        payload: t.Dict[str, t.Any] = py_json.loads(request.body or b"{}")
    except Exception as exc:  # noqa: BLE001 - keep simple for mocked endpoint
        return http_JsonResponse(data={"error": "invalid_json", "detail": str(exc)}, status=400)

    dataset = payload.get("dataset") or ""
    filters = payload.get("filters") or {}

    # Simple deterministic mock: base on dataset name hash and number of filters
    base = abs(hash(dataset)) % 10000
    multiplier = max(1, len(filters) or 1)
    mocked_count = 50000 + (base // 10) + (multiplier * 123)

    return http_JsonResponse(data={"count": mocked_count})


@http_decorators.require_GET
def cellinfo_filters_suggest(request):
    """
    Return mocked suggestions for a given string field and query.

    :param request: The HTTP request with query parameters: field, dataset, q, limit, filters.

    :raise: None

    :return: A JSON response containing a list of suggestion strings and a truncated flag.
    """
    field = request.GET.get("field", "")
    q = (request.GET.get("q", "") or "").strip()
    try:
        limit = int(request.GET.get("limit", "20"))
    except ValueError:
        limit = 20
    limit = max(1, min(limit, 50))

    # Very simple mock dictionary
    mock_vocab = {
        "cell_type": [
            "T cell",
            "T regulatory cell",
            "B cell",
            "Neuron",
            "Astrocyte",
            "Macrophage",
            "Endothelial cell",
        ],
        "assay": [
            "10x",
            "smart-seq",
            "Drop-seq",
            "microwell-seq",
            "BD Rhapsody",
        ],
    }

    candidates = mock_vocab.get(field, [])
    if q:
        q_lower = q.lower()
        candidates = [c for c in candidates if c.lower().startswith(q_lower)]

    suggestions = candidates[:limit]
    truncated = len(candidates) > len(suggestions)

    return http_JsonResponse(data={"suggestions": suggestions, "truncated": truncated})
