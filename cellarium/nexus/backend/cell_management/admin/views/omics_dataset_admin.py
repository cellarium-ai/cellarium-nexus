"""
Admin module for omics dataset management.
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
from cellarium.nexus.backend.cell_management.models import OmicsDataset, OmicsDatasetBackend
from cellarium.nexus.omics_datastore.bq_ops import create_bq_tables

logger = logging.getLogger(__name__)


@admin.register(OmicsDataset)
class OmicsDatasetAdmin(ModelAdmin):
    """
    Admin interface for managing omics datasets.

    Provides functionality to create and manage omics datasets with different backends.
    """

    list_display = ("id", "name", "backend", "schema", "description", "link_display")
    search_fields = ("name",)
    list_filter = ("name", "backend", "schema")
    ordering = ("name",)
    readonly_fields = ("link",)
    fieldsets = ((None, {"fields": ("name", "schema", "backend", "description", "uri", "link")}),)

    def save_model(self, request: HttpRequest, obj: OmicsDataset, form: forms.Form, change: bool) -> None:
        """
        Create backend resources when a new record is created.

        For BigQuery backend, create the BigQuery dataset in GCP.

        :param request: The HTTP request
        :param obj: The OmicsDataset instance being saved
        :param form: The form used to create/edit the instance
        :param change: Boolean indicating if this is a change to an existing record

        :raise Exception: If dataset creation fails
        """
        is_new = not change  # `change` is False if adding a new object

        if is_new and obj.backend == OmicsDatasetBackend.BIGQUERY:
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

        if is_new and obj.backend == OmicsDatasetBackend.BIGQUERY:
            html_link = format_html(constants.BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT, url_link=obj.link)
            self.message_user(
                request=request,
                message=format_html(constants.BIGQUERY_SUCCESS_MESSAGE_TEXT + " " + html_link),
                level=messages.SUCCESS,
            )

    def link_display(self, obj: OmicsDataset) -> str:
        """
        Generate clickable link for the dataset dashboard.

        :param obj: The current OmicsDataset object

        :return: HTML anchor tag with the link or plain text
        """
        if obj.link:
            return format_html(constants.BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT, url_link=obj.link)
        return "-"

    link_display.short_description = "Dashboard Link"
    link_display.allow_tags = True
