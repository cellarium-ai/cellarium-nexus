from django import forms
from django.utils.translation import gettext_lazy as _
from unfold.widgets import UnfoldAdminFileFieldWidget, UnfoldAdminSelectWidget, UnfoldAdminTextInputWidget
import json

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

    bigquery_dataset = forms.ModelChoiceField(
        label=_("BigQuery Dataset"),
        queryset=BigQueryDataset.objects.all(),
        widget=UnfoldAdminSelectWidget(attrs={'disabled': 'disabled'}),
        help_text=_("BigQuery Dataset to extract data from"),
        required=False,  # Not required since we handle it in the view
    )

    filters = forms.JSONField(
        label=_("Filters"),
        required=False,
        widget=forms.Textarea(attrs={'rows': 4}),
        help_text=_("JSON formatted filters to apply during extraction"),
    )

    def __init__(self, *args, **kwargs):
        """
        Initialize the form.

        If there's only one BigQuery dataset, make the field read-only.
        """
        super().__init__(*args, **kwargs)
        
        # If there's an initial dataset, make the field read-only
        if 'initial' in kwargs and kwargs['initial'].get('bigquery_dataset'):
            self.fields['bigquery_dataset'].widget.attrs['disabled'] = 'disabled'
            # Store the initial value to use it in clean
            self._initial_dataset = kwargs['initial'].get('bigquery_dataset')

    def clean(self):
        """
        Clean the form data.

        Ensure the BigQuery dataset is set, either from the form or from the initial value.
        """
        cleaned_data = super().clean()
        
        # If the field is disabled, it won't be in cleaned_data, so use the initial value
        if 'bigquery_dataset' not in cleaned_data and hasattr(self, '_initial_dataset'):
            cleaned_data['bigquery_dataset'] = self._initial_dataset
        
        return cleaned_data

    def clean_filters(self):
        """
        Validate the filters JSON field.

        :raise forms.ValidationError: If JSON is invalid

        :return: Parsed JSON filters or empty dict
        """
        filters = self.cleaned_data.get('filters')
        if not filters:
            return {}
        
        try:
            if isinstance(filters, str):
                return json.loads(filters)
            return filters
        except json.JSONDecodeError as e:
            raise forms.ValidationError(_("Invalid JSON format: %s") % str(e))
