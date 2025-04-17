from django.contrib.auth import get_user_model
from django.db import models
from django.utils.translation import gettext_lazy as _


UserModel = get_user_model()


def default_empty_dict():
    return {}


class Curriculum(models.Model):
    creator = models.ForeignKey(
        to=UserModel,
        on_delete=models.CASCADE,
        related_name="curriculums",
        verbose_name=_("creator"),
    )
    cell_count = models.IntegerField(verbose_name=_("cell count"))
    extract_bin_size = models.IntegerField(verbose_name=_("extract bin size"))
    extract_files_dir = models.CharField(
        max_length=512,
        verbose_name=_("extract files directory"),
    )
    metadata_file_path = models.CharField(
        max_length=512,
        verbose_name=_("metadata file path"),
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

    class Meta:
        verbose_name = _("curriculum")
        verbose_name_plural = _("curriculums")
        app_label = "curriculum"

    def __str__(self):
        return f"Curriculum {self.id} by {self.creator.username}"
