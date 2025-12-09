"""
Admin module for cell information management.
"""

import builtins
import decimal as py_decimal
import json as py_json
import logging
import typing as t

import django.forms as dj_forms
import django.urls as django_urls
import django.utils.html as django_html
import django.views.decorators.http as http_decorators
import django.views.generic as generic
import google.api_core as google_api_core
from django.conf import settings
from django.contrib import admin, messages
from django.http import HttpResponseRedirect as http_HttpResponseRedirect
from django.http import JsonResponse as http_JsonResponse
from google.cloud import bigquery

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.admin import forms as admin_forms
from cellarium.nexus.backend.cell_management.admin import schemas as admin_schemas
from cellarium.nexus.backend.cell_management.utils import bigquery_utils, soma_workflows_utils, workflows_utils
from cellarium.nexus.backend.cell_management.utils import filters as filters_utils
from cellarium.nexus.omics_datastore.bq_ops.bq_avro_schemas import cell_management as bq_schemas
from cellarium.nexus.omics_datastore.bq_ops.data_operator import BigQueryDataOperator
from cellarium.nexus.omics_datastore.protocols import DataOperatorProtocol
from cellarium.nexus.omics_datastore.soma_ops.data_operator import TileDBSOMADataOperator

logger = logging.getLogger(__name__)


def _create_cached_manager_for_dataset(*, dataset_name: str) -> bigquery_utils.OmicsCachedDataManager:
    """
    Create an OmicsCachedDataManager for a dataset based on its backend type.

    :param dataset_name: Name of the omics dataset

    :raise cell_models.OmicsDataset.DoesNotExist: If dataset not found

    :return: Configured OmicsCachedDataManager instance
    """
    dataset = cell_models.OmicsDataset.objects.get(name=dataset_name)
    operator: DataOperatorProtocol

    if dataset.backend == cell_models.OmicsDatasetBackend.TILEDB_SOMA:
        if not dataset.uri:
            raise ValueError(f"TileDB SOMA dataset '{dataset_name}' has no URI configured")
        operator = TileDBSOMADataOperator(experiment_uri=dataset.uri)
        cache_namespace = f"soma|{dataset.uri}"
    else:
        # Default to BigQuery
        bq_client = bigquery.Client(project=settings.GCP_PROJECT_ID)
        operator = BigQueryDataOperator(
            client=bq_client,
            project=settings.GCP_PROJECT_ID,
            dataset=dataset_name,
        )
        cache_namespace = f"bq|{settings.GCP_PROJECT_ID}|{dataset_name}"

    return bigquery_utils.OmicsCachedDataManager(
        operator=operator,
        cache_namespace=cache_namespace,
    )


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
        bq_datasets_qs = cell_models.OmicsDataset.objects.all().order_by("name")
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
        try:
            for ds in datasets:
                manager = _create_cached_manager_for_dataset(dataset_name=ds)
                dataset_counts[ds] = manager.get_cached_count(filters_dict={})
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
                    ds_manager = _create_cached_manager_for_dataset(dataset_name=ds)
                    cat = ds_manager.get_cached_categorical_obs_columns(
                        distinct_threshold=settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
                    )
                    precomputed_categorical_all[ds] = sorted(list(cat))
                    ds_map: dict[str, list[str]] = {}

                    # For SOMA, use all categorical columns; for BQ, filter through schema
                    dataset_obj = cell_models.OmicsDataset.objects.get(name=ds)
                    if dataset_obj.backend == cell_models.OmicsDatasetBackend.TILEDB_SOMA:
                        columns_to_fetch = cat
                    else:
                        columns_to_fetch = [col for col in string_columns if col in cat]

                    for col in columns_to_fetch:
                        vals = ds_manager.get_cached_distinct_obs_values(
                            column_name=col,
                        )
                        if vals:
                            ds_map[col] = vals
                    precomputed_suggestions_all[ds] = ds_map
        except Exception as exc:
            logger.exception(f"Failed to precompute categorical columns/suggestions for datasets {datasets}: {exc}")
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

    For BigQuery datasets, fields are derived from CellInfoBQAvroSchema.
    For TileDB SOMA datasets, fields are derived from categorical columns in the obs table.

    :param dataset: Optional dataset name to determine backend and categorical columns

    :return: List of field metadata dictionaries
    """
    results: list[dict] = []

    # Determine which string columns are categorical based on distinct count per dataset
    categorical_columns: set[str] = set()
    is_soma = False
    if dataset:
        try:
            dataset_obj = cell_models.OmicsDataset.objects.get(name=dataset)
            is_soma = dataset_obj.backend == cell_models.OmicsDatasetBackend.TILEDB_SOMA
            logger.info(
                f"_get_cellinfo_filters_fields: dataset={dataset}, backend={dataset_obj.backend}, is_soma={is_soma}"
            )
        except cell_models.OmicsDataset.DoesNotExist:
            logger.warning(f"_get_cellinfo_filters_fields: dataset={dataset} not found")

        manager = _create_cached_manager_for_dataset(dataset_name=dataset)
        categorical_columns = manager.get_cached_categorical_obs_columns(
            distinct_threshold=settings.FILTERS_CATEGORICAL_UNIQUE_LIMIT,
        )
        logger.info(f"_get_cellinfo_filters_fields: categorical_columns={categorical_columns}")

    if is_soma:
        # For SOMA, build fields directly from categorical columns (all are strings with suggestions)
        for col in sorted(categorical_columns):
            results.append(
                {
                    "key": col,
                    "label": col.replace("_", " ").title(),
                    "type": "string",
                    "operators": ["in", "not_in"],
                    "suggest": True,
                    "suggest_mode": "prefix",
                    "suggest_min_chars": 2,
                }
            )
    else:
        # For BigQuery, build fields from CellInfoBQAvroSchema
        schema_fields = bq_schemas.CellInfoBQAvroSchema.model_fields
        type_hints = t.get_type_hints(bq_schemas.CellInfoBQAvroSchema, include_extras=True)
        exclude_keys = {"metadata_extra"}

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
                    # suggest only for categorical string columns (<= limit distinct values)
                    suggest = key in categorical_columns if dataset else False
                    # Split operators for plain string vs categorical string
                    operators = ["in", "not_in"] if suggest else ["eq", "not_eq", "in", "not_in"]
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

    manager = _create_cached_manager_for_dataset(dataset_name=dataset)
    try:
        count = manager.get_cached_count(filters_dict=normalized_filters)
    except google_api_core.exceptions.GoogleAPICallError as exc:
        return http_JsonResponse(data={"error": "bigquery_error", "detail": str(exc)}, status=502)

    resp = admin_schemas.CountResponse(count=int(count))
    return http_JsonResponse(data=resp.model_dump())


class ExtractCurriculumAdminView(generic.FormView):
    """
    Render and submit the Extract Curriculum form.

    Pre-fill initial values from query parameters if provided:
    - ``dataset`` (omics dataset name)
    - ``filters`` (JSON-encoded string of filter statements)
    """

    template_name = "admin/cell_management/cellinfo/extract_curriculum.html"
    form_class = admin_forms.ExtractCurriculumForm

    def get_context_data(self, **kwargs):
        """
        Build template context and include admin site context for Unfold.

        :param kwargs: Keyword args from FormView

        :raise: None

        :return: Context dictionary for the template
        """
        context = super().get_context_data(**kwargs)
        context.update(admin.site.each_context(self.request))
        # Provide initial filters as a plain dict for template json_script
        try:
            form = context.get("form")
            initial_filters = {}
            if form is not None:
                initial_filters = form.initial.get("filters") or {}
            context["filters_initial"] = initial_filters
        except Exception:
            context["filters_initial"] = {}
        return context

    def get_initial(self) -> dict:
        """
        Build initial form data from query parameters.

        :raise: None

        :return: Initial dictionary for the form
        """
        initial: dict[str, t.Any] = super().get_initial()
        dataset_name = self.request.GET.get("dataset") or ""
        filters_raw = self.request.GET.get("filters") or ""

        if dataset_name:
            try:
                ds_obj = cell_models.OmicsDataset.objects.get(name=dataset_name)
                initial["omics_dataset"] = ds_obj.pk
            except cell_models.OmicsDataset.DoesNotExist:
                pass

        if filters_raw:
            try:
                parsed = py_json.loads(filters_raw)
            except Exception:
                parsed = {}
            # Normalize to backend format in case client passed loose structure
            try:
                normalized = filters_utils.normalize_filter_statements(filter_statements=parsed)
            except Exception:
                normalized = {}
            initial["filters"] = normalized

        return initial

    def form_valid(self, form: admin_forms.ExtractCurriculumForm):
        """
        Submit the extract pipeline using validated form data.

        Detect the dataset backend and submit the appropriate pipeline:
        - BigQuery datasets use the BQ extract pipeline
        - TileDB SOMA datasets use the SOMA extract pipeline

        :param form: Validated ExtractCurriculumForm

        :raise google.api_core.exceptions.GoogleAPIError: If pipeline submission fails
        :raise ValueError: If dataset configuration is invalid

        :return: HTTP redirect to the Vertex AI Pipeline run URL
        """
        cleaned = form.cleaned_data
        omics_dataset = cleaned["omics_dataset"]

        try:
            if omics_dataset.backend == cell_models.OmicsDatasetBackend.TILEDB_SOMA:
                # Submit SOMA extract pipeline
                pipeline_url = soma_workflows_utils.submit_soma_extract_pipeline(
                    name=cleaned["name"],
                    creator_id=self.request.user.id,
                    omics_dataset=omics_dataset,
                    range_size=settings.TILEDB_SOMA_EXTRACT_CONTIGUOUS_RANGE,
                    output_chunk_size=cleaned["extract_bin_size"],
                    filters=cleaned.get("filters") or None,
                    shuffle_ranges=True,
                    feature_schema=cleaned.get("feature_schema"),
                    obs_columns=cleaned.get("obs_columns") or None,
                    var_columns=settings.TILEDB_SOMA_EXTRACT_VAR_COLUMNS,
                    x_layer=settings.TILEDB_SOMA_EXTRACT_X_LAYER,
                    output_format=settings.TILEDB_SOMA_EXTRACT_OUTPUT_FORMAT,
                    max_workers_extract=16,
                    max_workers_shuffle=6,
                )
            else:
                # Submit BigQuery extract pipeline
                pipeline_url = workflows_utils.submit_extract_pipeline(
                    feature_schema=cleaned["feature_schema"],
                    creator_id=self.request.user.id,
                    omics_dataset=omics_dataset,
                    name=cleaned["name"],
                    extract_bin_size=cleaned["extract_bin_size"],
                    categorical_column_count_limit=cleaned["categorical_column_count_limit"],
                    obs_columns=cleaned["obs_columns"],
                    extract_bin_keys=cleaned.get("extract_bin_keys") or None,
                    filters=cleaned.get("filters") or None,
                    metadata_extra_columns=cleaned.get("metadata_extra_columns") or None,
                )
            messages.success(
                self.request,
                django_html.format_html(
                    'Extract pipeline submitted successfully. <a href="{url}" target="_blank" rel="noopener noreferrer" style="text-decoration: underline;">View run</a>',
                    url=pipeline_url,
                ),
            )
        except Exception as exc:  # Surface error to admin as message and redisplay form
            messages.error(self.request, f"Failed to run extract pipeline: {exc}")
            return self.form_invalid(form)

        # Redirect back to Cell Info page within Admin
        return http_HttpResponseRedirect(django_urls.reverse("admin:cell_management_cellinfo_changelist"))
