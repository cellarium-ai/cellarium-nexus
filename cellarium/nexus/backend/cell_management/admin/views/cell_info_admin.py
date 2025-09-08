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
from django.conf import settings
import builtins
import decimal as py_decimal
import json as py_json
import typing as t

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.admin.utils import bigquery_utils
from cellarium.nexus.backend.cell_management.admin import forms as admin_forms
from cellarium.nexus.backend.cell_management.admin import schemas as admin_schemas
from cellarium.nexus.backend.cell_management.admin.views.utils import filters as filters_utils
from cellarium.nexus.omics_datastore.bq_avro_schemas import cell_management as bq_schemas
from cellarium.nexus.omics_datastore.bq_ops import constants as bq_constants


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
        datasets_list: list[dict[str, str]] = [
            {"name": ds.name, "description": ds.description or ""} for ds in bq_datasets_qs
        ]
        datasets: list[str] = [ds["name"] for ds in datasets_list]

        # Choose a default selection: use singleton default if present, else first
        default_ds_obj = (
            bq_datasets_qs.get_default_dataset() if hasattr(bq_datasets_qs, "get_default_dataset") else None
        )
        selected_dataset = default_ds_obj.name if default_ds_obj else (datasets[0] if datasets else "")

        # Compute counts with caching; handle BigQuery errors gracefully
        dataset_counts: dict[str, int] = {}
        manager = bigquery_utils.BigQueryCachedDataManager()
        try:
            for ds in datasets:
                dataset_counts[ds] = manager.get_cached_count_bq(dataset_name=ds, filters_dict={})
        except google_api_core.exceptions.GoogleAPICallError:
            messages.error(
                request=self.request,
                message=(
                    "Unable to retrieve cell counts from BigQuery at the moment. "
                    "Please try again later or contact an administrator."
                ),
            )

        # Build server-rendered filters formset
        fields_meta = _get_cellinfo_filters_fields(dataset=selected_dataset)
        field_choices: list[tuple[str, str]] = [(f["key"], f.get("label") or f["key"]) for f in fields_meta]

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

        # Precompute across all datasets for instant switching
        precomputed_categorical_all: dict[str, list[str]] = {}
        precomputed_suggestions_all: dict[str, dict[str, list[str]]] = {}
        try:
            if selected_dataset:
                # Determine string columns from schema
                schema_fields = bq_schemas.CellInfoBQAvroSchema.model_fields
                type_hints = t.get_type_hints(bq_schemas.CellInfoBQAvroSchema, include_extras=True)
                string_columns: list[str] = []
                for key, field_info in schema_fields.items():
                    if key == "metadata_extra" or key.endswith("_ontology_term_id"):
                        continue
                    anno = type_hints.get(key, field_info.annotation)
                    base = filters_utils.resolve_base_type(anno)
                    if base in (str,):
                        string_columns.append(key)

                for ds in datasets:
                    if not ds:
                        continue
                    cat = manager.get_cached_categorical_columns_bq(
                        dataset_name=ds,
                        table_name=bq_constants.BQ_CELL_INFO_TABLE_NAME,
                        distinct_threshold=settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
                    )
                    precomputed_categorical_all[ds] = sorted(list(cat))
                    ds_map: dict[str, list[str]] = {}
                    for col in string_columns:
                        if col not in cat:
                            continue
                        vals = manager.get_cached_distinct_values_bq(
                            dataset_name=ds,
                            column_name=col,
                            table_name=bq_constants.BQ_CELL_INFO_TABLE_NAME,
                            limit=settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
                        )
                        if vals:
                            ds_map[col] = vals
                    precomputed_suggestions_all[ds] = ds_map
        except Exception:
            # Do not block page rendering if suggestions warm-up fails
            precomputed_categorical_all = {}
            precomputed_suggestions_all = {}

        context.update(
            {
                "datasets": datasets_list,
                "selected_dataset": selected_dataset,
                "dataset_counts": dataset_counts,
                "filters_formset": formset,
                "filters_fields_meta": fields_meta,
                "precomputed_categorical_columns_all": precomputed_categorical_all,
                "precomputed_suggestions_all": precomputed_suggestions_all,
                **admin.site.each_context(self.request),
            }
        )
        return context


def _get_cellinfo_filters_fields(*, dataset: str | None = None) -> list[dict]:
    """
    Return filters metadata used by the server-rendered formset and API endpoints.

    :raise: None

    :return: List of field metadata dictionaries
    """
    # Dynamically build fields from CellInfoBQAvroSchema, excluding metadata_extra and ontology term ids
    schema_fields = bq_schemas.CellInfoBQAvroSchema.model_fields
    type_hints = t.get_type_hints(bq_schemas.CellInfoBQAvroSchema, include_extras=True)
    exclude_keys = {"metadata_extra"}
    results: list[dict] = []

    # Determine which string columns are categorical based on distinct count per dataset
    # Use cached manager instead of querying directly from the view
    categorical_columns: set[str] = set()
    if dataset:
        manager = bigquery_utils.BigQueryCachedDataManager()
        categorical_columns = manager.get_cached_categorical_columns_bq(
            dataset_name=dataset,
            table_name=bq_constants.BQ_CELL_INFO_TABLE_NAME,
            distinct_threshold=settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
        )

    for key, field_info in schema_fields.items():
        if key in exclude_keys or key.endswith("_ontology_term_id"):
            continue

        # Determine type and operators
        anno = type_hints.get(key, field_info.annotation)
        ftype: str
        operators: list[str]
        suggest: bool = False

        base = filters_utils.resolve_base_type(anno)

        match base:
            # --- scalar primitives (note: bool is a subclass of int; handle it first)
            case builtins.bool:
                ftype = "boolean"
                operators = ["eq", "not_eq"]
            case builtins.float | py_decimal.Decimal | builtins.int:
                ftype = "number"
                operators = ["eq", "gt", "gte", "lt", "lte"]
            case builtins.str:
                ftype = "string"
                operators = ["eq", "not_eq", "in", "not_in"]
                # suggest only for categorical string columns (<= limit distinct values)
                suggest = key in categorical_columns if dataset else False
            case _:
                # Skip unsupported/complex types
                continue

        results.append(
            {
                "key": key,
                "label": (field_info.title or key.replace("_", " ").title()),
                "type": ftype,
                "operators": operators,
                "suggest": suggest,
                **({"suggest_mode": "prefix", "suggest_min_chars": 2} if suggest else {}),
            }
        )

    return results


@http_decorators.require_POST
def cellinfo_filters_count(request):
    """
    Return a mocked count based on provided dataset and filters.

    :param request: The HTTP request.

    :raise ValueError: If payload is not valid JSON.

    :return: A JSON response containing a mocked count integer.
    """
    try:
        raw: t.Dict[str, t.Any] = py_json.loads(request.body or b"{}")
        parsed = admin_schemas.FiltersPayload(**raw)
    except Exception as exc:  # noqa: BLE001 - keep simple for mocked endpoint
        return http_JsonResponse(data={"error": "invalid_json", "detail": str(exc)}, status=400)

    dataset = (parsed.dataset or "").strip()
    if not dataset:
        return http_JsonResponse(data={"error": "invalid_dataset", "detail": "Dataset is required."}, status=400)

    filter_statements = parsed.filters
    normalized_filters = filters_utils.normalize_filter_statements(filter_statements=filter_statements)

    manager = bigquery_utils.BigQueryCachedDataManager()
    try:
        count = manager.get_cached_count_bq(dataset_name=dataset, filters_dict=normalized_filters)
    except google_api_core.exceptions.GoogleAPICallError as exc:
        return http_JsonResponse(data={"error": "bigquery_error", "detail": str(exc)}, status=502)

    resp = admin_schemas.CountResponse(count=int(count))
    return http_JsonResponse(data=resp.model_dump())
