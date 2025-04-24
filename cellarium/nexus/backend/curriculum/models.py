from django.contrib.auth import get_user_model
from django.db import models
from django.utils.translation import gettext_lazy as _

UserModel = get_user_model()


def default_empty_dict():
    return {}


class Curriculum(models.Model):
    STATUS_PREPARE = "PREPARE"
    STATUS_EXTRACTING = "EXTRACTING"
    STATUS_SUCCEEDED = "SUCCEEDED"
    STATUS_FAILED = "FAILED"
    STATUS_CHOICES = [
        (STATUS_PREPARE, "Prepare initiated..."),
        (STATUS_EXTRACTING, "Extracting..."),
        (STATUS_SUCCEEDED, "Succeeded"),
        (STATUS_FAILED, "Failed"),
    ]

    name = models.CharField(max_length=512, verbose_name=_("name"))
    creator = models.ForeignKey(
        to=UserModel,
        on_delete=models.CASCADE,
        related_name="curriculums",
        verbose_name=_("creator"),
    )
    cell_count = models.IntegerField(verbose_name=_("cell count"), null=True, blank=True)
    extract_bin_size = models.IntegerField(verbose_name=_("extract bin size"), null=True, blank=True)
    extract_bin_count = models.IntegerField(verbose_name=_("extract bin count"), null=True, blank=True)
    extract_files_path = models.FilePathField(
        verbose_name=_("extract files directory"),
        null=True,
        blank=True,
    )
    metadata_file_path = models.FileField(
        verbose_name=_("metadata file path"),
        null=True,
        blank=True,
    )
    created_at = models.DateTimeField(
        verbose_name=_("created at"),
        auto_now_add=True,
        editable=False,
    )
    filters_json = models.JSONField(
        verbose_name=_("filters json"),
        default=default_empty_dict,
        null=True,
        blank=True,
    )
    status = models.CharField(
        max_length=20,
        choices=STATUS_CHOICES,
        default=STATUS_PREPARE,
        verbose_name=_("status"),
    )

    class Meta:
        verbose_name = _("curriculum")
        verbose_name_plural = _("curriculums")
        app_label = "curriculum"

    def __str__(self):
        return f"Curriculum {self.id} by {self.creator.username}"
