from django.apps import AppConfig
from django.utils.translation import gettext_lazy as _


class IngestManagementConfig(AppConfig):
    verbose_name = _("Ingest Management")
    name = "cellarium.nexus.backend.ingest_management"
    label = "ingest_management"
