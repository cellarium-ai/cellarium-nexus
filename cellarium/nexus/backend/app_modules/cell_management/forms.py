import json

from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from nexus.omics_datastore.bq_avro_schemas.column_mapping import ColumnMappingSchema
from unfold.widgets import UnfoldAdminFileFieldWidget, UnfoldAdminSelectWidget, UnfoldAdminTextInputWidget

from .models import BigQueryDataset, ColumnMapping, FeatureSchema


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


class CreateSchemaFromCSVForm(forms.Form):
    """
    Form for creating a new schema from CSV file.
    """

    schema_name = forms.CharField(
        label=_("Schema Name"),
        max_length=255,
        help_text=_("Name for the new feature schema"),
        widget=UnfoldAdminTextInputWidget,
    )
    csv_file = forms.FileField(
        label=_("Features CSV"),
        help_text=_("CSV file with ensemble_id and symbol columns"),
        widget=UnfoldAdminFileFieldWidget,
    )


class PrepareExtractTablesForm(forms.Form):
    """
    Form for preparing extract tables.
    """

    feature_schema = forms.ModelChoiceField(
        label=_("Feature Schema"),
        queryset=FeatureSchema.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("Feature schema to use for extraction"),
    )

    extract_table_prefix = forms.CharField(
        label=_("Extract Table Prefix"),
        max_length=255,
        help_text=_("Prefix for the extract tables"),
        widget=UnfoldAdminTextInputWidget,
    )
