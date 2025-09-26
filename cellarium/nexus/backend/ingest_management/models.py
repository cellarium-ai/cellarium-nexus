import uuid

from django.contrib.auth import get_user_model
from django.db import models
from django.utils.translation import gettext_lazy as _

from cellarium.nexus.backend.ingest_management.utils.column_mapping_utils import (
    get_obs_column_choices,
    get_var_column_choices,
)


class IngestInfo(models.Model):
    STATUS_STARTED = "STARTED"
    STATUS_SUCCEEDED = "SUCCEEDED"
    STATUS_FAILED = "FAILED"
    STATUS_CHOICES = [(STATUS_STARTED, "Started"), (STATUS_SUCCEEDED, "Succeeded"), (STATUS_FAILED, "Failed")]
    nexus_uuid = models.UUIDField(default=uuid.uuid4, verbose_name=_("nexus uuid"), unique=True)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default=STATUS_STARTED, verbose_name=_("status"))
    bigquery_dataset = models.ForeignKey(
        to="cell_management.BigQueryDataset",
        on_delete=models.CASCADE,
        related_name="ingest_management_ingests",
        verbose_name=_("BigQuery dataset"),
    )
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
        verbose_name = _("ingest info")
        verbose_name_plural = _("ingest info objects")
        app_label = "ingest_management"

    def __str__(self):
        return f"Ingest {str(self.id)} info"


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
    bigquery_dataset = models.ForeignKey(
        to="cell_management.BigQueryDataset",
        on_delete=models.CASCADE,
        related_name="ingest_management_index_trackings",
        verbose_name=_("BigQuery dataset"),
    )
    resource_key = models.CharField(max_length=100, verbose_name=_("resource key"))
    largest_index = models.BigIntegerField(default=0, verbose_name=_("largest assigned index"))

    class Meta:
        verbose_name = "Index Tracking"
        verbose_name_plural = "Index Tracking Entries"
        app_label = "ingest_management"
        unique_together = ("bigquery_dataset", "resource_key")

    def __str__(self):
        return f"{self.bigquery_dataset} | {self.resource_key} -- {self.largest_index}"


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
    input_file_gcs_path = models.CharField(max_length=1024, verbose_name=_("input file GCS path"))
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

        # Truncate the GCS path if it's too long
        max_length = 40
        gcs_path = self.input_file_gcs_path
        if gcs_path and len(gcs_path) > max_length:
            # Extract the filename (last part of the path)
            path_parts = gcs_path.split("/")
            filename = path_parts[-1] if path_parts else ""

            # Use the filename or a truncated version with ellipsis
            if len(filename) <= max_length:
                truncated_path = filename
            else:
                truncated_path = f"{gcs_path[:15]}...{gcs_path[-20:]}"
        else:
            truncated_path = gcs_path

        return f"{status} - {self.validator_name} - {truncated_path}"
