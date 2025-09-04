"""
Admin module for cell information management.
"""

from django.contrib import admin
from django.contrib import messages
import django.views.generic as generic
import google.api_core as google_api_core

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.cell_management.admin.utils import bigquery_utils


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

        context.update({
            "datasets": datasets,
            "selected_dataset": selected_dataset,
            "dataset_counts": dataset_counts,
            **admin.site.each_context(self.request),
        })
        return context
