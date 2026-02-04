import uuid

from django.contrib.auth import get_user_model
from django.db import models
from django.utils.translation import gettext_lazy as _

from cellarium.nexus.backend.ingest_management.utils.column_mapping_utils import (
    get_obs_column_choices,
    get_var_column_choices,
)

STATUS_STARTED = "STARTED"
STATUS_SUCCEEDED = "SUCCEEDED"
STATUS_FAILED = "FAILED"
STATUS_CHOICES = [
    (STATUS_STARTED, "Started"),
    (STATUS_SUCCEEDED, "Succeeded"),
    (STATUS_FAILED, "Failed"),
]


def default_empty_list():
    return []


class Ingest(models.Model):
    STATUS_STARTED = STATUS_STARTED
    STATUS_SUCCEEDED = STATUS_SUCCEEDED
    STATUS_FAILED = STATUS_FAILED
    STATUS_CHOICES = STATUS_CHOICES

    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default=STATUS_STARTED, verbose_name=_("status"))
    omics_dataset = models.ForeignKey(
        to="cell_management.OmicsDataset",
        on_delete=models.CASCADE,
        related_name="ingest_management_ingests",
        verbose_name=_("omics dataset"),
    )
    metadata_extra = models.JSONField(verbose_name=_("metadata extra"), null=True, blank=True)
    ingest_start_timestamp = models.DateTimeField(
        verbose_name=_("ingest start timestamp"), auto_now_add=True, editable=False
    )
    ingest_finish_timestamp = models.DateTimeField(verbose_name=_("ingest finish timestamp"), null=True, blank=True)
    gencode_version = models.PositiveSmallIntegerField(verbose_name=_("gencode version"), null=True, blank=True)

    class Meta:
        verbose_name = _("ingest")
        verbose_name_plural = _("ingests")
        ordering = ["-ingest_start_timestamp"]
        app_label = "ingest_management"

    def __str__(self):
        return f"Ingest {str(self.id)}"


class IngestInfo(models.Model):
    STATUS_STARTED = STATUS_STARTED
    STATUS_SUCCEEDED = STATUS_SUCCEEDED
    STATUS_FAILED = STATUS_FAILED
    STATUS_CHOICES = STATUS_CHOICES

    ingest = models.ForeignKey(
        to="ingest_management.Ingest",
        on_delete=models.CASCADE,
        related_name="files",
        verbose_name=_("ingest"),
    )
    nexus_uuid = models.UUIDField(default=uuid.uuid4, verbose_name=_("nexus uuid"), unique=True)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default=STATUS_STARTED, verbose_name=_("status"))
    omics_dataset = models.ForeignKey(
        to="cell_management.OmicsDataset",
        on_delete=models.CASCADE,
        related_name="ingest_management_ingest_files",
        verbose_name=_("omics dataset"),
    )
    gcs_file_path = models.CharField(
        max_length=1024,
        verbose_name=_("gcs file path"),
        null=True,
        blank=True,
    )
    tag = models.CharField(max_length=255, verbose_name=_("tag"), null=True, blank=True)
    metadata_extra = models.JSONField(verbose_name=_("metadata extra"), null=True, blank=True)
    ingest_start_timestamp = models.DateTimeField(
        verbose_name=_("ingest start timestamp"), auto_now_add=True, editable=False
    )
    ingest_finish_timestamp = models.DateTimeField(verbose_name=_("ingest finish timestamp"), null=True, blank=True)
    column_mapping = models.JSONField(
        verbose_name=_("column mapping"),
        null=True,
        blank=True,
        help_text=_("Mapping of input column names to schema column names"),
    )

    class Meta:
        verbose_name = _("ingest file")
        verbose_name_plural = _("ingest files")
        app_label = "ingest_management"

    def __str__(self):
        return f"Ingest file {str(self.id)}"


class ColumnMapping(models.Model):
    """
    Model for storing reusable column mappings for data ingestion.
    """

    name = models.CharField(max_length=255, verbose_name=_("name"), unique=True)
    description = models.TextField(verbose_name=_("description"), null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True, verbose_name=_("created at"))
    updated_at = models.DateTimeField(auto_now=True, verbose_name=_("updated at"))

    class Meta:
        verbose_name = _("column mapping")
        verbose_name_plural = _("column mappings")
        ordering = ["-updated_at"]
        app_label = "ingest_management"

    def __str__(self):
        return self.name


class ObsColumnMapping(models.Model):
    """
    Model for storing obs column mappings.
    """

    mapping = models.ForeignKey(ColumnMapping, on_delete=models.CASCADE, related_name="obs_mappings")
    input_column = models.CharField(max_length=255, verbose_name=_("input column"))
    schema_column = models.CharField(max_length=255, verbose_name=_("schema column"), choices=get_obs_column_choices)

    class Meta:
        verbose_name = _("obs column mapping")
        verbose_name_plural = _("obs column mappings")
        unique_together = [("mapping", "input_column"), ("mapping", "schema_column")]
        app_label = "ingest_management"

    def __str__(self):
        return f"{self.input_column} -> {self.schema_column}"


class VarColumnMapping(models.Model):
    """
    Model for storing var column mappings.
    """

    mapping = models.ForeignKey(ColumnMapping, on_delete=models.CASCADE, related_name="var_mappings")
    input_column = models.CharField(max_length=255, verbose_name=_("input column"))
    schema_column = models.CharField(max_length=255, verbose_name=_("schema column"), choices=get_var_column_choices)

    class Meta:
        verbose_name = _("var column mapping")
        verbose_name_plural = _("var column mappings")
        unique_together = [("mapping", "input_column"), ("mapping", "schema_column")]
        app_label = "ingest_management"

    def __str__(self):
        return f"{self.input_column} -> {self.schema_column}"


class IndexTracking(models.Model):
    omics_dataset = models.ForeignKey(
        to="cell_management.OmicsDataset",
        on_delete=models.CASCADE,
        related_name="ingest_management_index_trackings",
        verbose_name=_("omics dataset"),
    )
    resource_key = models.CharField(max_length=100, verbose_name=_("resource key"))
    largest_index = models.BigIntegerField(default=0, verbose_name=_("largest assigned index"))

    class Meta:
        verbose_name = "Index Tracking"
        verbose_name_plural = "Index Tracking Entries"
        app_label = "ingest_management"
        unique_together = ("omics_dataset", "resource_key")

    def __str__(self):
        return f"{self.omics_dataset} | {self.resource_key} -- {self.largest_index}"


class ValidationReport(models.Model):
    """
    Model for storing validation reports for data ingestion.
    """

    created_at = models.DateTimeField(auto_now_add=True, verbose_name=_("created at"))
    creator = models.ForeignKey(
        to=get_user_model(),
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="validation_reports",
        verbose_name=_("creator"),
    )

    class Meta:
        verbose_name = _("validation report")
        verbose_name_plural = _("validation reports")
        ordering = ["-created_at"]
        app_label = "ingest_management"

    def __str__(self):
        return f"Validation Report {self.id} - {self.created_at}"


class ValidationReportItem(models.Model):
    """
    Model for storing individual validation report items.
    """

    report = models.ForeignKey(
        ValidationReport,
        on_delete=models.CASCADE,
        related_name="items",
        verbose_name=_("report"),
    )
    input_file_path = models.CharField(max_length=1024, verbose_name=_("input file path"))
    sanitized_file_path = models.CharField(
        max_length=1024,
        verbose_name=_("sanitized file path"),
        null=True,
        blank=True,
    )
    validator_name = models.CharField(max_length=255, verbose_name=_("validator name"))
    is_valid = models.BooleanField(verbose_name=_("is valid"))
    message = models.TextField(verbose_name=_("message"), null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True, verbose_name=_("created at"))

    class Meta:
        verbose_name = _("validation report item")
        verbose_name_plural = _("validation report items")
        ordering = ["-created_at"]
        app_label = "ingest_management"

    def __str__(self):
        """
        Return a string representation of the validation report item with truncated GCS path.

        :return: String representation
        """
        status = "Valid" if self.is_valid else "Invalid"

        # Truncate the path if it's too long
        max_length = 40
        file_path = self.input_file_path
        if file_path and len(file_path) > max_length:
            # Extract the filename (last part of the path)
            path_parts = file_path.split("/")
            filename = path_parts[-1] if path_parts else ""

            # Use the filename or a truncated version with ellipsis
            if len(filename) <= max_length:
                truncated_path = filename
            else:
                truncated_path = f"{file_path[:15]}...{file_path[-20:]}"
        else:
            truncated_path = file_path

        return f"{status} - {self.validator_name} - {truncated_path}"


class IngestSchema(models.Model):
    """
    Model for storing ingest schemas.
    """

    class XValidationType(models.TextChoices):
        COUNT_MATRIX = "count_matrix", _("Count matrix")
        FEATURE_MATRIX = "feature_matrix", _("Feature matrix")

    name = models.CharField(max_length=255, verbose_name=_("name"), unique=True)
    description = models.TextField(verbose_name=_("description"), null=True, blank=True)
    x_validation_type = models.CharField(
        max_length=32,
        choices=XValidationType.choices,
        default=XValidationType.COUNT_MATRIX,
        verbose_name=_("X validation type"),
    )
    measurement_name = models.CharField(
        max_length=256,
        verbose_name=_("measurement name"),
        help_text=_("Type of measurement (e.g., RNA, ATAC, OPS)"),
    )
    created_at = models.DateTimeField(auto_now_add=True, verbose_name=_("created at"))
    updated_at = models.DateTimeField(auto_now=True, verbose_name=_("updated at"))

    class Meta:
        verbose_name = _("ingest schema")
        verbose_name_plural = _("ingest schemas")
        ordering = ["-updated_at"]
        app_label = "ingest_management"

    def __str__(self):
        return self.name


class SomaVarSchema(models.Model):
    """
    Model for storing SOMA var schema definitions as Parquet files.
    """

    ingest_schema = models.OneToOneField(
        "IngestSchema",
        on_delete=models.CASCADE,
        related_name="var_schema",
        verbose_name=_("ingest schema"),
    )
    is_subset = models.BooleanField(
        default=True,
        verbose_name=_("is subset"),
        help_text=_("If true, input AnnData may contain a subset of features from the schema"),
    )
    var_parquet_file = models.FileField(
        upload_to="soma/var_schemas/",
        verbose_name=_("var parquet file"),
    )
    feature_count = models.PositiveIntegerField(
        default=0,
        verbose_name=_("feature count"),
        help_text=_("Number of features in the schema (derived from parquet index)"),
    )
    var_columns = models.JSONField(
        default=default_empty_list,
        verbose_name=_("var columns"),
        help_text=_("List of var metadata columns stored in the parquet file"),
    )
    created_at = models.DateTimeField(auto_now_add=True, verbose_name=_("created at"))
    updated_at = models.DateTimeField(auto_now=True, verbose_name=_("updated at"))

    class Meta:
        verbose_name = _("SOMA var schema")
        verbose_name_plural = _("SOMA var schemas")
        ordering = ["-updated_at"]
        app_label = "ingest_management"

    def __str__(self):
        return f"{self.ingest_schema}"


class SomaObsColumnSchema(models.Model):
    """
    Model for storing SOMA obs column schema entries.
    """

    class DType(models.TextChoices):
        BOOL = "bool", _("bool")
        INT8 = "int8", _("int8")
        INT16 = "int16", _("int16")
        INT32 = "int32", _("int32")
        INT64 = "int64", _("int64")
        FLOAT32 = "float32", _("float32")
        FLOAT64 = "float64", _("float64")
        STRING = "string", _("string")
        CATEGORY = "category", _("category")

    ingest_schema = models.ForeignKey(
        IngestSchema,
        on_delete=models.CASCADE,
        related_name="obs_columns",
        verbose_name=_("ingest schema"),
    )
    name = models.CharField(max_length=255, verbose_name=_("name"))
    dtype = models.CharField(max_length=32, choices=DType.choices, verbose_name=_("dtype"))
    nullable = models.BooleanField(default=False, verbose_name=_("nullable"))

    class Meta:
        verbose_name = _("SOMA obs column schema")
        verbose_name_plural = _("SOMA obs column schemas")
        app_label = "ingest_management"
        unique_together = [("ingest_schema", "name")]

    def __str__(self):
        return f"{self.name} ({self.dtype})"
