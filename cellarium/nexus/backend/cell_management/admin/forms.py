import json

from django import forms
from django.forms.widgets import TextInput
from django.utils.translation import gettext_lazy as _
from django_json_widget.widgets import JSONEditorWidget as BaseJSONEditorWidget
from unfold.widgets import (
    UnfoldAdminFileFieldWidget,
    UnfoldAdminIntegerFieldWidget,
    UnfoldAdminSelectWidget,
    UnfoldAdminTextInputWidget,
)

from cellarium.nexus.backend.cell_management.models import BigQueryDataset, FeatureSchema
from cellarium.nexus.backend.curriculum.models import Curriculum


class CommaSeparatedWidget(UnfoldAdminTextInputWidget):
    """
    Widget that displays a list of strings as comma-separated values.
    """

    def __init__(self, attrs=None):
        """
        Initialize the widget with placeholder text.

        :param attrs: HTML attributes
        """
        default_attrs = {
            "class": "vTextField",
            "placeholder": "Enter column names separated by commas",
            "style": "width: 100%;",
        }
        if attrs:
            default_attrs.update(attrs)
        super().__init__(attrs=default_attrs)

    def format_value(self, value):
        """
        Format the value for display in the widget.

        :param value: Value to format

        :return: Formatted value
        """
        if value is None:
            return ""
        if isinstance(value, list):
            return ", ".join(str(v) for v in value)
        return value


class CommaSeparatedField(forms.CharField):
    """
    Field that converts comma-separated input to a list of strings.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the field with the CommaSeparatedWidget.

        :param args: Positional arguments
        :param kwargs: Keyword arguments
        """
        kwargs.setdefault("widget", CommaSeparatedWidget)
        super().__init__(*args, **kwargs)

    def to_python(self, value):
        """
        Convert the input value to a list of strings.

        :param value: Input value

        :raise: forms.ValidationError

        :return: List of strings
        """
        if not value:
            return []
        if isinstance(value, list):
            return value
        return [item.strip() for item in value.split(",") if item.strip()]


class CustomJSONEditorWidget(BaseJSONEditorWidget):
    """
    Custom JSON editor widget that properly handles both string and dict inputs.
    """

    def __init__(
        self,
        attrs: dict | None = None,
        mode: str = "tree",
        options: dict | None = None,
        width: str | None = None,
        height: str | None = None,
        **kwargs,
    ):
        """
        Initialize the widget.

        :param attrs: HTML attributes for the widget
        :param mode: Editor mode (tree, code, etc.)
        :param options: Additional options for the editor
        :param width: Editor width
        :param height: Editor height
        :param kwargs: Additional keyword arguments
        """
        if options is None:
            options = {}
        if attrs is None:
            attrs = {}

        if width:
            attrs["style"] = f"width: {width}"
        if height:
            attrs["style"] = attrs.get("style", "") + f"; height: {height}"

        super().__init__(attrs=attrs, mode=mode, options=options, **kwargs)

    def format_value(self, value):
        """
        Format the value for display in the widget.

        :param value: Value to format
        :return: Formatted value
        """
        if isinstance(value, str):
            try:
                return json.loads(value)
            except json.JSONDecodeError:
                return value
        return value


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


class ExtractCurriculumForm(forms.Form):
    """
    Form for submitting extract curriculum job.
    """

    feature_schema = forms.ModelChoiceField(
        label=_("Feature Schema"),
        queryset=FeatureSchema.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("Feature schema to use for extraction"),
    )
    name = forms.CharField(
        label=_("Curriculum Name"),
        max_length=255,
        help_text=_("Name for the new curriculum"),
        widget=UnfoldAdminTextInputWidget,
    )
    extract_bin_size = forms.IntegerField(
        label=_("Extract Bin Size"),
        min_value=1,
        max_value=100_000,
        help_text=_("Bin size for the extract tables"),
        widget=UnfoldAdminIntegerFieldWidget,
    )
    bigquery_dataset = forms.ModelChoiceField(
        label=_("BigQuery Dataset"),
        queryset=BigQueryDataset.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("BigQuery Dataset to extract data from"),
    )
    metadata_extra_columns = CommaSeparatedField(
        label=_("Extra Metadata Columns"),
        required=False,
        help_text=_("Enter column names separated by commas to include as additional metadata in the extract."),
    )
    filters = forms.JSONField(
        label=_("Filters"),
        required=False,
        widget=CustomJSONEditorWidget(
            attrs={
                "height": "400px",
                "width": "100%",
            },
            options={
                "modes": ["tree", "code"],
                "mode": "tree",
                "search": True,
                "sortObjectKeys": True,
                "enableSort": False,
                "enableTransform": False,
            },
        ),
        help_text=_("JSON formatted filters to apply during extraction. Use tree view for easier editing."),
    )

    def __init__(self, *args, **kwargs):
        """
        Initialize the form.

        If there's only one BigQuery dataset, make the field read-only.
        """
        # If initial data contains filters as a string, parse it to dict
        if "initial" in kwargs and "filters" in kwargs["initial"]:
            filters = kwargs["initial"]["filters"]
            if isinstance(filters, str):
                try:
                    kwargs["initial"]["filters"] = json.loads(filters)
                except json.JSONDecodeError:
                    # If parsing fails, keep as is
                    pass

        super().__init__(*args, **kwargs)

    def clean_name(self):
        """
        Validate that no curriculum with this name exists.

        :raise: forms.ValidationError

        :return: Validated name
        """
        name = self.cleaned_data.get("name")
        if Curriculum.objects.filter(name=name).exists():
            raise forms.ValidationError(_("A curriculum with this name already exists"))
        return name

    def clean_filters(self):
        """
        Validate the filters JSON field.

        :raise: forms.ValidationError

        :return: Parsed JSON filters or empty dict
        """
        filters = self.cleaned_data.get("filters")
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

    def clean_metadata_extra_columns(self):
        """
        Validate the metadata_extra_columns field.

        :return: List of column name strings
        """
        columns = self.cleaned_data.get("metadata_extra_columns")
        return columns if columns else []
