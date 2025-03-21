from django.apps import AppConfig
from django.utils.translation import gettext_lazy as _


class CurriculumConfig(AppConfig):
    """
    Configuration for the curriculum app.
    """
    default_auto_field = "django.db.models.BigAutoField"
    name = "nexus.backend.app_modules.curriculum"
    verbose_name = _("Curriculum") 