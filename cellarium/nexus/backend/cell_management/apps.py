from django.apps import AppConfig
from django.utils.translation import gettext_lazy as _


class CellManagementConfig(AppConfig):
    name = "cellarium.nexus.backend.cell_management"
    label = "cell_management"
    default_auto_field = "django.db.models.BigAutoField"
    verbose_name = _("Cell Management")
