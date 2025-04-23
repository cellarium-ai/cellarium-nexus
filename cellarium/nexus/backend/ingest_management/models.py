import uuid

from django.contrib.contenttypes.fields import GenericForeignKey
from django.db import models
from django.utils.translation import gettext_lazy as _

from cellarium.nexus.backend.ingest_management.utils.column_mapping import (
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
    content_type = models.OneToOneField(
        to="contenttypes.ContentType",
        on_delete=models.CASCADE,
        unique=True,
        verbose_name=_("tracked model"),
        related_name="ingest_management_indextracking",
    )
    table_object_id = models.PositiveIntegerField(default=0, verbose_name=_("table object id"))
    tracked_model = GenericForeignKey(ct_field="content_type", fk_field="table_object_id")
    largest_index = models.BigIntegerField(default=0, verbose_name=_("largest assigned index"))

    class Meta:
        verbose_name = "Index Tracking"
        verbose_name_plural = "Index Tracking Entries"
        app_label = "ingest_management"

    def __str__(self):
        return f"{self.content_type.name} -- {self.largest_index}"
