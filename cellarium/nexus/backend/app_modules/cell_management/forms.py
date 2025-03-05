from django import forms
from django.utils.translation import gettext_lazy as _
from unfold.widgets import UnfoldAdminFileFieldWidget, UnfoldAdminSelectWidget, UnfoldAdminTextInputWidget
from django_json_widget.widgets import JSONEditorWidget as BaseJSONEditorWidget
import json

from .models import BigQueryDataset, ColumnMapping, FeatureSchema


class CustomJSONEditorWidget(BaseJSONEditorWidget):
    """
    Custom JSON editor widget that properly handles both string and dict inputs.
    """
    def __init__(self, attrs: dict | None = None, mode: str = 'tree', options: dict | None = None, width: str | None = None, height: str | None = None) -> None:
        default_options = {
            'modes': ['tree', 'code'],
            'mode': mode,
            'search': True,
            'onLoad': '''function (editor) {
                editor.expandAll();
            }'''
        }
        if options:
            default_options.update(options)

        super().__init__(attrs=attrs, mode=mode, options=default_options, width=width, height=height)

    def format_value(self, value: dict | str | None) -> dict | None:
        """
        Format the value for the widget.

        :param value: The value to format (can be string or dict)
        :return: Python dict or None
        """
        if value is None:
            return None
            
        if isinstance(value, dict):
            return value
            
        if isinstance(value, str):
            try:
                return json.loads(value)
            except json.JSONDecodeError:
                return None
                
        return None


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
        widget=CustomJSONEditorWidget(
            attrs={
                'height': '400px',
                'width': '100%',
            },
            options={
                'modes': ['tree', 'code'],
                'mode': 'tree',
                'search': True,
                'sortObjectKeys': True,
                'enableSort': False,  
                'enableTransform': False,
            }
        ),
        help_text=_("JSON formatted filters to apply during extraction. Use tree view for easier editing."),
    )

    def __init__(self, *args, **kwargs):
        """
        Initialize the form.

        If there's only one BigQuery dataset, make the field read-only.
        """
        # If initial data contains filters as a string, parse it to dict
        if 'initial' in kwargs and 'filters' in kwargs['initial']:
            filters = kwargs['initial']['filters']
            if isinstance(filters, str):
                try:
                    kwargs['initial']['filters'] = json.loads(filters)
                except json.JSONDecodeError:
                    # If parsing fails, keep as is
                    pass
        
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
        
        # If filters is already a dict, return it
        if isinstance(filters, dict):
            return filters
            
        try:
            if isinstance(filters, str):
                return json.loads(filters)
            return filters
        except json.JSONDecodeError as e:
            raise forms.ValidationError(_("Invalid JSON format: %s") % str(e))
