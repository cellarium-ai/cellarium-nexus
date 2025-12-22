from django.db import models
from django.utils.translation import gettext_lazy as _


def default_empty_dict():
    return {}


class OmicsDatasetBackend(models.TextChoices):
    """Backend types for omics datasets."""

    BIGQUERY = "bigquery", _("BigQuery")
    TILEDB_SOMA = "tiledb_soma", _("TileDB SOMA")


class OmicsDatasetQuerySet(models.QuerySet):
    def get_default_dataset(self) -> models.Model | None:
        """
        Get the default omics dataset if only one exists.

        :return: The default omics dataset or None if none or multiple exist
        """
        dataset_count = self.model.objects.count()

        if dataset_count == 1:
            dataset = self.model.objects.first()
            if dataset:
                return dataset

        return None


class OmicsDataset(models.Model):
    """Model for storing omics dataset configurations."""

    name = models.CharField(max_length=256, verbose_name=_("name"), unique=True)
    description = models.TextField(verbose_name=_("description"), null=True, blank=True)
    backend = models.CharField(
        max_length=32,
        choices=OmicsDatasetBackend.choices,
        default=OmicsDatasetBackend.BIGQUERY,
        verbose_name=_("backend"),
        help_text=_("The storage backend for this dataset"),
    )
    link = models.CharField(
        max_length=512,
        verbose_name=_("link"),
        unique=True,
        null=True,
        blank=True,
        help_text=_("Link to the Google Cloud dashboard for this dataset"),
    )
    uri = models.CharField(
        max_length=512,
        verbose_name=_("URI"),
        unique=True,
        null=True,
        blank=True,
        help_text=_("URI for the dataset (e.g., gs:// path for SOMA, or BigQuery dataset name)"),
    )
    objects = OmicsDatasetQuerySet.as_manager()

    class Meta:
        verbose_name = _("omics dataset")
        verbose_name_plural = _("omics datasets")
        app_label = "cell_management"

    def __str__(self):
        return self.name


# Backwards compatibility alias for migrations
BigQueryDataset = OmicsDataset


class FeatureSchema(models.Model):
    """
    Model for storing feature schemas and their associated features.
    """

    name = models.CharField(max_length=255, verbose_name=_("name"), unique=True)

    class Meta:
        verbose_name = _("feature schema")
        verbose_name_plural = _("feature schemas")
        app_label = "cell_management"

    def __str__(self):
        return self.name


class FeatureInfo(models.Model):
    """
    Model for storing feature information like ensemble IDs and symbols.

    Each feature belongs to exactly one schema and will be deleted when the schema is deleted.
    Each ensemble_id must be unique within a schema.
    """

    ensemble_id = models.CharField(max_length=255, verbose_name=_("ensemble id"))
    symbol = models.CharField(max_length=255, verbose_name=_("symbol"), blank=True, null=True)
    feature_schema = models.ForeignKey(
        to="FeatureSchema",
        on_delete=models.CASCADE,
        related_name="features",
        verbose_name=_("feature schema"),
    )

    class Meta:
        verbose_name = _("feature info")
        verbose_name_plural = _("feature info objects")
        app_label = "cell_management"
        unique_together = [("feature_schema", "ensemble_id")]
        indexes = [
            models.Index(fields=["feature_schema", "ensemble_id"]),
        ]

    def __str__(self):
        return f"{self.symbol} ({self.ensemble_id})"
