from django import forms
from django.utils.translation import gettext_lazy as _
from unfold.widgets import UnfoldAdminFileFieldWidget, UnfoldAdminSelectWidget

from cellarium.nexus.backend.cell_management.models import BigQueryDataset
from cellarium.nexus.backend.ingest_management.models import ColumnMapping


class ValidateNewDataChangeListActionForm(forms.Form):
    """
    Form for validating new data uploads.
    """

    ingest_csv_file = forms.FileField(
        label=_("Input datasets"),
        widget=UnfoldAdminFileFieldWidget,
        help_text=_("CSV file with GCS Bucket paths of files to ingest"),
    )


class IngestNewDataChangeListActionForm(forms.Form):
    """
    Form for ingesting new data.
    """

    ingest_csv_file = forms.FileField(
        label=_("Input datasets"),
        widget=UnfoldAdminFileFieldWidget,
        help_text=_("CSV file with GCS Bucket paths of files to ingest"),
    )
    bigquery_dataset = forms.ModelChoiceField(
        label=_("BigQuery Dataset"),
        queryset=BigQueryDataset.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("BigQuery Dataset to ingest the input files to"),
    )
    column_mapping = forms.ModelChoiceField(
        required=False,
        label=_("Column Mapping"),
        queryset=ColumnMapping.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("Select an existing column mapping configuration"),
    )

    def clean_column_mapping(self):
        """
        Validate the column mapping model instance if provided.

        :return: ColumnMapping instance or None
        """
        mapping = self.cleaned_data.get("column_mapping")
        if not mapping:
            return None

        return mapping
