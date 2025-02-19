from django.apps import AppConfig
from django.utils.translation import gettext_lazy as _


class CellManagementConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    verbose_name = _("Cell Management")
    name = "nexus.backend.app_modules.cell_management"

    def ready(self):
        from nexus.backend.app_modules.cell_management.api import views  # noqa
