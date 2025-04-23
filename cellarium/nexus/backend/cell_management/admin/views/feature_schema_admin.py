"""
Admin module for feature schema management.
"""

import csv
import logging

import pandas as pd
from django.contrib import admin, messages
from django.db.models import QuerySet
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import redirect, render
from django.urls import path, reverse
from django.utils.html import format_html, format_html_join
from django.utils.translation import gettext_lazy as _
from unfold.admin import ModelAdmin
from unfold.decorators import action

from cellarium.nexus.backend.cell_management.admin import constants, forms
from cellarium.nexus.backend.cell_management.models import FeatureInfo, FeatureSchema

logger = logging.getLogger(__name__)


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
        form = forms.CreateSchemaFromCSVForm(data=request.POST or None, files=request.FILES or None)

        if request.method == "POST" and form.is_valid():
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

        return render(
            request=request,
            template_name=constants.CHANGELIST_ACTION_FORM,
            context={
                "form": form,
                "title": _("Upload Features CSV"),
                "submit_button_title": _("Upload"),
                **self.admin_site.each_context(request),
            },
        )
