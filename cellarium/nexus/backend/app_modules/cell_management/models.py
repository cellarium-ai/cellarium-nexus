import uuid

from django.contrib.contenttypes.fields import GenericForeignKey
from django.db import models
from django.utils.translation import gettext_lazy as _
from nexus.backend.app_modules.cell_management.utils.column_mapping import (
    get_obs_column_choices,
    get_var_column_choices,
)


def default_empty_dict():
    return {}


class BigQueryDataset(models.Model):
    name = models.CharField(max_length=256, verbose_name=_("name"), unique=True)
    description = models.TextField(verbose_name=_("description"), null=True, blank=True)
    link = models.CharField(max_length=512, verbose_name=_("link"), unique=True, null=True, blank=True)

    class Meta:
        verbose_name = _("BigQuery dataset")
        verbose_name_plural = _("BigQuery datasets")
        app_label = "cell_management"

    def __str__(self):
        return self.name


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
        related_name="ingests",
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

    def __str__(self):
        return f"Ingest {str(self.id)} info"


class CellFeatureInfo(models.Model):
    id = models.BigAutoField(primary_key=True, verbose_name=_("ID"), editable=False)
    ensemble_id = models.CharField(max_length=256, verbose_name=_("ensemble id"))
    symbol = models.CharField(max_length=255, verbose_name=_("symbol"))
    biotype = models.CharField(max_length=255, blank=True, null=True, verbose_name=_("biotype"))
    is_filtered = models.BooleanField(blank=True, null=True, verbose_name=_("is filtered"))
    reference = models.CharField(max_length=255, blank=True, null=True, verbose_name=_("reference"))
    ingest = models.ForeignKey(
        to="cell_management.IngestInfo", on_delete=models.CASCADE, related_name="features", verbose_name=_("ingest")
    )
    tag = models.CharField(max_length=64, verbose_name=_("tag"), null=True, db_index=True)
    metadata_extra = models.JSONField(
        verbose_name=_("metadata extra"), default=default_empty_dict, null=True, blank=True
    )

    class Meta:
        verbose_name = _("cell feature info")
        verbose_name_plural = _("cell feature info objects")
        indexes = [models.Index(fields=["ensemble_id"])]

    def __str__(self):
        return f"{self.symbol} ({self.ensemble_id})"


class CellInfo(models.Model):
    id = models.BigAutoField(primary_key=True, verbose_name=_("ID"), editable=False)
    original_id = models.CharField(max_length=256, verbose_name=_("original id"), unique=True)
    ingest = models.ForeignKey(
        to="cell_management.IngestInfo", on_delete=models.CASCADE, related_name="cells", verbose_name=_("ingest")
    )
    tag = models.CharField(max_length=64, verbose_name=_("tag"), null=True, db_index=True)
    metadata_extra = models.JSONField(
        verbose_name=_("metadata extra"), default=default_empty_dict, null=True, blank=True
    )
    # Cell Features
    donor_id = models.CharField(max_length=256, null=True, blank=True, verbose_name=_("donor id"))
    cell_type = models.CharField(max_length=256, null=False, db_index=True, verbose_name=_("cell type"))
    assay = models.CharField(max_length=256, null=True, blank=True, db_index=True, verbose_name=_("assay"))
    development_stage = models.CharField(
        max_length=256, null=True, blank=True, db_index=True, verbose_name=_("development stage")
    )
    tissue = models.CharField(max_length=256, null=True, blank=True, db_index=True, verbose_name=_("tissue"))
    disease = models.CharField(max_length=256, null=True, blank=True, db_index=True, verbose_name=_("disease"))
    organism = models.CharField(max_length=256, null=True, blank=True, db_index=True, verbose_name=_("organism"))
    self_reported_ethnicity = models.CharField(
        max_length=256, null=True, blank=True, db_index=True, verbose_name=_("self-reported ethnicity")
    )
    sex = models.CharField(max_length=256, null=True, blank=True, db_index=True, verbose_name=_("sex"))
    suspension_type = models.CharField(
        max_length=256, null=True, blank=True, db_index=True, verbose_name=_("suspension type")
    )
    total_mrna_umis = models.IntegerField(null=True, blank=True, db_index=True, verbose_name=_("total mrna umis"))

    # Cell Features Ontology Term IDs
    cell_type_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, db_index=True, verbose_name=_("cell type ontology term id")
    )
    assay_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("assay ontology term id")
    )
    development_stage_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("development stage ontology term id")
    )
    tissue_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("tissue ontology term id")
    )
    disease_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("disease ontology term id")
    )
    organism_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("organism ontology term id")
    )
    self_reported_ethnicity_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("self-reported ethnicity ontology term id")
    )
    sex_ontology_term_id = models.CharField(
        max_length=256, null=True, blank=True, verbose_name=_("sex ontology term id")
    )

    class Meta:
        verbose_name = _("cell info")
        verbose_name_plural = _("cel info objects")

    def __str__(self):
        return f"Cell {str(self.id)[:6]} info"


class IndexTracking(models.Model):
    content_type = models.OneToOneField(
        to="contenttypes.ContentType", on_delete=models.CASCADE, unique=True, verbose_name=_("tracked model")
    )
    table_object_id = models.PositiveIntegerField(default=0, verbose_name=_("table object id"))
    tracked_model = GenericForeignKey(ct_field="content_type", fk_field="table_object_id")
    largest_index = models.BigIntegerField(default=0, verbose_name=_("largest assigned index"))

    class Meta:
        verbose_name = "Index Tracking"
        verbose_name_plural = "Index Tracking Entries"

    def __str__(self):
        return f"{self.content_type.name} -- {self.largest_index}"


class FeatureInfo(models.Model):
    """
    Model for storing feature information like ensemble IDs and symbols.
    """

    ensemble_id = models.CharField(max_length=255, verbose_name=_("ensemble id"), unique=True)
    symbol = models.CharField(max_length=255, verbose_name=_("symbol"))

    class Meta:
        verbose_name = _("feature info")
        verbose_name_plural = _("feature info objects")
        app_label = "cell_management"

    def __str__(self):
        return f"{self.symbol} ({self.ensemble_id})"


class FeatureSchema(models.Model):
    """
    Model for storing feature schemas and their associated features.
    """

    name = models.CharField(max_length=255, verbose_name=_("name"), unique=True)
    features = models.ManyToManyField(
        to="cell_management.FeatureInfo", related_name="schemas", verbose_name=_("features")
    )

    class Meta:
        verbose_name = _("feature schema")
        verbose_name_plural = _("feature schemas")

    def __str__(self):
        return self.name


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

    def __str__(self) -> str:
        return f"{self.input_column} → {self.schema_column}"


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

    def __str__(self) -> str:
        return f"{self.input_column} → {self.schema_column}"
