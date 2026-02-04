from django import forms
from django.utils.translation import gettext_lazy as _
from unfold.widgets import (
    UnfoldAdminFileFieldWidget,
    UnfoldAdminSelect2Widget,
    UnfoldAdminSelectWidget,
    UnfoldAdminTextInputWidget,
)

from cellarium.nexus.backend.cell_management.models import OmicsDataset
from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.utils import soma_schema_utils

GENCODE_CHOICES = [
    (43, _("Gencode 43")),
    (44, _("Gencode 44")),
]

GENCODE_CHOICES_WITH_NONE = [(None, _("---------"))] + GENCODE_CHOICES


class IngestNewDataChangeListActionForm(forms.Form):
    """
    Form for ingesting new data.
    """

    ingest_csv_file = forms.FileField(
        label=_("Input datasets"),
        widget=UnfoldAdminFileFieldWidget,
        help_text=_("CSV file with GCS Bucket paths of files to ingest"),
    )
    omics_dataset = forms.ModelChoiceField(
        label=_("Omics Dataset"),
        queryset=OmicsDataset.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("Omics dataset to ingest the input files to"),
    )
    column_mapping = forms.ModelChoiceField(
        required=False,
        label=_("Column Mapping"),
        queryset=models.ColumnMapping.objects.all(),
        widget=UnfoldAdminSelectWidget,
        help_text=_("Select an existing column mapping configuration"),
    )
    gencode_version = forms.ChoiceField(
        required=False,
        initial=None,
        label=_("Gencode Version"),
        choices=GENCODE_CHOICES_WITH_NONE,
        widget=UnfoldAdminSelectWidget,
        help_text=_(
            "Select the Gencode version for validation. Will be ignored for non Homo Sapiens data. If not provided, "
            "will not be used for validation."
        ),
    )

    def clean_column_mapping(self) -> models.ColumnMapping | None:
        """
        Validate the column mapping model instance if provided.

        :return: ColumnMapping instance or None
        """
        mapping = self.cleaned_data.get("column_mapping")
        if not mapping:
            return None

        return mapping

    def clean_gencode_version(self) -> int | None:
        gencode_version = self.cleaned_data.get("gencode_version")
        if gencode_version == "":
            return None

        return gencode_version


class ValidateNewDataChangeListActionForm(forms.Form):
    """
    Form for SOMA validation and sanitization pipeline.
    """

    omics_dataset = forms.ModelChoiceField(
        label=_("Omics Dataset"),
        queryset=OmicsDataset.objects.all(),
        widget=UnfoldAdminSelect2Widget,
        help_text=_("Omics dataset containing the schema for validation"),
    )
    input_csv_file = forms.FileField(
        label=_("Input H5AD Files CSV"),
        widget=UnfoldAdminFileFieldWidget,
        help_text=_("CSV file with one GCS path per line (e.g., gs://bucket/file1.h5ad)"),
    )
    output_directory_uri = forms.CharField(
        label=_("Output Directory URI"),
        widget=UnfoldAdminTextInputWidget,
        help_text=_("GCS directory URI for sanitized output files (e.g., gs://bucket/output/)"),
        required=True,
    )

    def clean(self):
        cleaned = super().clean()
        output_directory = cleaned.get("output_directory_uri", "").strip()

        if output_directory and not output_directory.startswith("gs://"):
            self.add_error("output_directory_uri", _("Output directory URI must start with gs://"))

        return cleaned


class SomaIngestForm(forms.Form):
    """
    Form for SOMA ingest from validated and sanitized files.
    """

    omics_dataset = forms.ModelChoiceField(
        label=_("Omics Dataset"),
        queryset=OmicsDataset.objects.all(),
        widget=UnfoldAdminSelect2Widget,
        help_text=_("Omics dataset to ingest into"),
    )
    ingest_batch_size = forms.IntegerField(
        label=_("Ingest Batch Size"),
        initial=10,
        min_value=1,
        widget=UnfoldAdminTextInputWidget,
        help_text=_("Number of h5ad files per partition (for parallel processing)"),
    )
    measurement_name = forms.CharField(
        label=_("Measurement Name"),
        initial="RNA",
        max_length=256,
        widget=UnfoldAdminTextInputWidget,
        help_text=_("Name of the SOMA measurement (e.g., RNA, ATAC)"),
    )


class SomaVarSchemaInlineForm(forms.ModelForm):
    """
    Inline form for creating a SOMA var schema from a CSV file.
    """

    csv_file = forms.FileField(
        label=_("Var Schema CSV"),
        required=True,
        widget=UnfoldAdminFileFieldWidget,
        help_text=_("CSV where the first column will be used as the index."),
    )

    class Meta:
        model = models.SomaVarSchema
        fields = ("is_subset", "csv_file")
        widgets = {}

    def clean(self):
        cleaned = super().clean()
        csv_file = cleaned.get("csv_file")
        if not csv_file:
            return cleaned
        if csv_file:
            try:
                parsed = soma_schema_utils.parse_var_schema_csv(csv_file=csv_file)
            except soma_schema_utils.VarSchemaCSVError as exc:
                self.add_error("csv_file", str(exc))
                return cleaned

            self._parsed_df = parsed.dataframe
            self._parsed_from_csv = True

        return cleaned
