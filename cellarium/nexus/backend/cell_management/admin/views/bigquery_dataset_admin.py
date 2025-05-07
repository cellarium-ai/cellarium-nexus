"""
Admin module for BigQuery dataset management.
"""

import logging

from django import forms
from django.conf import settings
from django.contrib import admin, messages
from django.http import HttpRequest
from django.utils.html import format_html
from google.cloud import bigquery
from unfold.admin import ModelAdmin

from cellarium.nexus.backend.cell_management.admin import constants
from cellarium.nexus.backend.cell_management.models import BigQueryDataset
from cellarium.nexus.backend.core.admin.helpers import CountRelatedObjectsDeleteMixin
from cellarium.nexus.omics_datastore.bq_ops import create_bq_tables

logger = logging.getLogger(__name__)


@admin.register(BigQueryDataset)
class BigQueryDatasetAdmin(CountRelatedObjectsDeleteMixin, ModelAdmin):
    """
    Admin interface for managing BigQuery datasets.

    Provides functionality to create and manage BigQuery datasets in GCP.
    """

    # Override the actions to use our custom delete_selected action
    actions = ["delete_selected"]

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
                bq_client = bigquery.Client()
                link_to_dataset = create_bq_tables.create_bigquery_objects(
                    client=bq_client,
                    project=settings.GCP_PROJECT_ID,
                    dataset=obj.name,
                    location="US",
                    labels={"application": settings.GCP_APPLICATION_BILLING_LABEL},
                )
                obj.link = link_to_dataset
            except Exception as e:
                self.message_user(request, f"Failed to create BigQuery dataset: {e}", level=messages.ERROR)
                return

        super().save_model(request=request, obj=obj, form=form, change=change)

        if is_new:
            html_link = format_html(constants.BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT, url_link=obj.link)
            self.message_user(
                request=request,
                message=format_html(constants.BIGQUERY_SUCCESS_MESSAGE_TEXT + " " + html_link),
                level=messages.SUCCESS,
            )

    def link_display(self, obj: BigQueryDataset) -> str:
        """
        Generate clickable link to BigQuery dataset.

        :param obj: The current BigQueryDataset object

        :return: HTML anchor tag with the link
        """
        if obj.link:
            return format_html(constants.BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT, url_link=obj.link)
        return "-"

    link_display.short_description = "BigQuery Link"
    link_display.allow_tags = True
