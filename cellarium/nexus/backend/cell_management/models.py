from django.db import models
from django.utils.translation import gettext_lazy as _


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


class CellFeatureInfo(models.Model):
    id = models.BigAutoField(primary_key=True, verbose_name=_("ID"), editable=False)
    ensemble_id = models.CharField(max_length=256, verbose_name=_("ensemble id"))
    symbol = models.CharField(max_length=255, verbose_name=_("symbol"))
    biotype = models.CharField(max_length=255, blank=True, null=True, verbose_name=_("biotype"))
    is_filtered = models.BooleanField(blank=True, null=True, verbose_name=_("is filtered"))
    reference = models.CharField(max_length=255, blank=True, null=True, verbose_name=_("reference"))
    ingest = models.ForeignKey(
        to="ingest_management.IngestInfo", on_delete=models.CASCADE, related_name="features", verbose_name=_("ingest")
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
        to="ingest_management.IngestInfo", on_delete=models.CASCADE, related_name="cells", verbose_name=_("ingest")
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
