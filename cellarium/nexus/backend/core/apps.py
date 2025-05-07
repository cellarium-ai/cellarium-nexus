from django.apps import AppConfig
from django.utils.translation import gettext_lazy as _


class CoreConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "cellarium.nexus.backend.core"
    label = "cellarium_backend_core"
    verbose_name = _("Core")

    def ready(self):
        from django.conf import settings  # noqa

        from cellarium.nexus.workflows import kubeflow  # noqa

        kubeflow.set_configs(base_image=settings.PIPELINE_BASE_IMAGE, service_account=settings.PIPELINE_SERVICE_ACCOUNT)
