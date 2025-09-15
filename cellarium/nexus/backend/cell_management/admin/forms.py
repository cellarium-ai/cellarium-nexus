import json

from django import forms
from django.utils.html import escape
from django.utils.translation import gettext_lazy as _
from django_json_widget.widgets import JSONEditorWidget as BaseJSONEditorWidget
from unfold import widgets as unfold_widgets

from cellarium.nexus.backend.cell_management.admin import constants
from cellarium.nexus.backend.cell_management.admin.utils import check_curriculum_exists
from cellarium.nexus.backend.cell_management.models import BigQueryDataset, FeatureSchema
from cellarium.nexus.backend.curriculum.models import Curriculum


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
        widget=unfold_widgets.UnfoldAdminTextInputWidget,
    )
    csv_file = forms.FileField(
        label=_("Features CSV"),
        help_text=_("CSV file with ensemble_id and symbol columns"),
        widget=unfold_widgets.UnfoldAdminFileFieldWidget,
    )


class UnfoldReadOnlyJSONWidget(unfold_widgets.UnfoldAdminTextareaWidget):
    """
    Theme-aware (Unfold) read-only JSON:
      • rows = number of JSON lines (clamped)
      • scrolls inside the field
      • selectable/copyable
    """

    min_rows = 3
    max_rows = 40

    def render(self, name, value, attrs=None, renderer=None):
        attrs = attrs or {}

        # Parse + pretty print
        if value in (None, ""):
            value = {}
        if isinstance(value, str):
            try:
                value = json.loads(value)
            except Exception:
                # it's a raw string, show as-is
                pretty = value
            else:
                pretty = json.dumps(value, indent=2, ensure_ascii=False)
        else:
            pretty = json.dumps(value, indent=2, ensure_ascii=False)

        # rows = number of lines (clamped)
        line_count = pretty.count("\n") + 1
        rows = max(self.min_rows, min(self.max_rows, line_count))
        attrs["rows"] = rows
        attrs["readonly"] = True

        style = attrs.get("style", "")
        style += ";resize:none;overflow:auto;max-height:28rem"
        attrs["style"] = style.lstrip(";")

        return super().render(name, escape(pretty), attrs, renderer)


class FreeListField(forms.MultipleChoiceField):
    def validate(self, value):
        # skip default check against self.choices
        pass


class ExtractCurriculumForm(forms.Form):
    """
    Form for submitting extract curriculum job.
    """

    feature_schema = forms.ModelChoiceField(
        label=_("Feature Schema"),
        queryset=FeatureSchema.objects.all(),
        widget=unfold_widgets.UnfoldAdminSelect2Widget,
        help_text=_("Feature schema to use for extraction"),
    )
    name = forms.CharField(
        label=_("Curriculum Name"),
        max_length=255,
        help_text=_("Name for the new curriculum"),
        widget=unfold_widgets.UnfoldAdminTextInputWidget,
    )
    extract_bin_size = forms.IntegerField(
        label=_("Extract Bin Size"),
        min_value=1,
        max_value=100_000,
        initial=10_000,
        help_text=_("Bin size for the extract tables"),
        widget=unfold_widgets.UnfoldAdminIntegerFieldWidget,
    )
    bigquery_dataset = forms.ModelChoiceField(
        label=_("BigQuery Dataset"),
        queryset=BigQueryDataset.objects.all(),
        widget=unfold_widgets.UnfoldAdminSelect2Widget,
        help_text=_("BigQuery Dataset to extract data from"),
    )
    obs_columns = forms.MultipleChoiceField(
        label=_("Observation Columns"),
        choices=constants.CELL_INFO_EXTRACT_COLUMNS_CHOICES,
        widget=unfold_widgets.UnfoldAdminSelect2MultipleWidget(
            attrs={
                "data-placeholder": "Select columns…",
                "data-minimum-results-for-search": "0",
                "data-allow-clear": "true",
                "data-close-on-select": "false",
            }
        ),
    )

    metadata_extra_columns = FreeListField(
        label=_("Extra Metadata Columns"),
        required=False,
        help_text=_("Enter column names separated by commas to include as additional metadata in the extract."),
        widget=unfold_widgets.UnfoldAdminSelect2MultipleWidget(
            attrs={
                "data-tags": "true",
                "data-placeholder": "Start typing…",
                "data-close-on-select": "false",
                "data-allow-clear": "true",
                "data-minimum-results-for-search": "0",
                "data-token-separators": ",",
            }
        ),
    )
    categorical_column_count_limit = forms.IntegerField(
        label=_("Categorical Column Count Limit"),
        min_value=1,
        max_value=25_000,
        initial=2000,
        widget=unfold_widgets.UnfoldAdminIntegerFieldWidget,
        help_text=_(
            "Maximum number of categories per categorical column to be considered as categorical. "
            "If the number of categories exceeds this limit, the column will not be unified across all extract files."
        ),
    )
    extract_bin_keys = forms.MultipleChoiceField(
        label=_("Extract Bin Keys"),
        required=False,
        help_text=_(
            "Optional list of comma-separated keys to bin by. If not provided, bins will be assigned randomly."
        ),
        choices=constants.CELL_INFO_EXTRACT_BIN_KEYS_COLUMNS_CHOICES,
        widget=unfold_widgets.UnfoldAdminSelect2MultipleWidget(
            attrs={
                "data-placeholder": "Select columns…",
                "data-minimum-results-for-search": "0",
                "data-allow-clear": "true",
                "data-close-on-select": "false",
            }
        ),
    )
    filters = forms.Field(
        required=False,
        label=_("Filters"),
        help_text=_(
            "JSON-formatted filters used during extraction. Read-only. To make changes, go back and choose "
            "the proper filter set."
        ),
        widget=UnfoldReadOnlyJSONWidget(),
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
        Validate that no curriculum with this name exists in the database or in the GCS bucket.

        :raise: forms.ValidationError: If a curriculum with this name already exists
        :raise: google.cloud.exceptions.GoogleCloudError: If there's an error communicating with Google Cloud Storage

        :return: Validated name
        """
        name = self.cleaned_data.get("name")

        # Check if curriculum exists in the database
        if Curriculum.objects.filter(name=name).exists():
            raise forms.ValidationError(_("A curriculum with this name already exists in the database"))

        # Check if curriculum exists in the GCS bucket
        if check_curriculum_exists(name=name):
            raise forms.ValidationError(_("A curriculum with this name already exists in the storage bucket"))

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

    class Media:
        css = {"all": ("admin/cell_management/cellinfo/css/extract.css",)}
        js: tuple[str, ...] = ()


class FilterRowForm(forms.Form):
    """
    Render a single filter row for the Cell Info admin page.

    :param args: Positional arguments forwarded to ``forms.Form``
    :param kwargs: Keyword arguments forwarded to ``forms.Form``

    :raise: None

    :return: None
    """

    # Field choices and operator choices are provided at runtime through ``__init__``
    field = forms.ChoiceField(
        label=_("Field"),
        choices=(),
        widget=unfold_widgets.UnfoldAdminSelectWidget(
            attrs={
                "style": "width: 100%;",
                "data-theme": "admin-autocomplete",
                "class": "unfold-admin-autocomplete admin-autocomplete",
            }
        ),
        required=True,
    )
    operator = forms.ChoiceField(
        label=_("Operator"),
        choices=(),
        widget=unfold_widgets.UnfoldAdminSelectWidget(
            attrs={
                "style": "width: 100%;",
                "data-theme": "admin-autocomplete",
                "class": "unfold-admin-autocomplete admin-autocomplete",
            }
        ),
        required=True,
    )
    value = forms.CharField(
        label=_("Value"),
        widget=unfold_widgets.UnfoldAdminTextInputWidget(
            attrs={
                "style": "width: 100%;",
                "placeholder": "Enter value",
            }
        ),
        required=True,
    )

    # Prototype-only fields to expose proper Unfold widgets and Media for cloning in JS.
    # These are not validated (required=False) and are rendered hidden in the template.
    value_text = forms.CharField(
        label=_("Value (text)"),
        widget=unfold_widgets.UnfoldAdminTextInputWidget(
            attrs={
                "style": "width: 100%;",
                "class": "vTextField w-full h-9 px-3 py-2 rounded border border-base-200 dark:border-base-700 bg-white dark:bg-base-900 text-base-700 dark:text-base-200 placeholder-base-400 focus:outline-none focus:ring focus:ring-primary-300 focus:border-primary-600 shadow-sm appearance-none",
                "placeholder": "Enter value",
            }
        ),
        required=False,
    )

    value_categorical = forms.ChoiceField(
        label=_("Value (categorical)"),
        choices=(),
        widget=unfold_widgets.UnfoldAdminSelectWidget(
            attrs={
                "style": "width: 100%;",
                "data-theme": "admin-autocomplete",
                "class": "unfold-admin-autocomplete admin-autocomplete",
            }
        ),
        required=False,
    )

    value_boolean = forms.ChoiceField(
        label=_("Value (boolean)"),
        choices=(
            ("", "Select value"),
            ("true", "true"),
            ("false", "false"),
        ),
        widget=unfold_widgets.UnfoldAdminSelectWidget(
            attrs={
                "style": "width: 100%;",
                "data-theme": "admin-autocomplete",
                "class": "unfold-admin-autocomplete admin-autocomplete",
            }
        ),
        required=False,
    )

    def __init__(self, *args, **kwargs):
        """
        Initialize the filter row form with dynamic choices.

        ``field`` choices come from ``kwargs['field_choices']`` if provided, otherwise remain empty.
        ``operator`` choices default to a union of supported operators unless ``operator_choices`` is provided.

        :param args: Positional arguments
        :param kwargs: Keyword arguments; may include ``field_choices`` and ``operator_choices``

        :raise: None

        :return: None
        """
        field_choices = kwargs.pop("field_choices", ())
        operator_choices = kwargs.pop("operator_choices", ())
        super().__init__(*args, **kwargs)

        if field_choices:
            self.fields["field"].choices = field_choices

        if operator_choices:
            self.fields["operator"].choices = operator_choices
        else:
            # Provide a sensible default union of operators used in the UI
            default_ops = [
                ("eq", _("equals")),
                ("not_eq", _("not equals")),
                ("in", _("in")),
                ("not_in", _("not in")),
                ("gt", ">"),
                ("gte", "\u2265"),
                ("lt", "<"),
                ("lte", "\u2264"),
            ]
            self.fields["operator"].choices = default_ops

    class Media:
        js = (
            "admin/js/vendor/jquery/jquery.js",
            "admin/js/vendor/select2/select2.full.js",
            "admin/js/jquery.init.js",
            "unfold/js/select2.init.js",
        )
        css = {
            "screen": (
                "admin/css/vendor/select2/select2.css",
                "admin/css/autocomplete.css",
            )
        }
