from django.contrib import admin
from django.urls import include, path
from django.views.generic import RedirectView

from cellarium.nexus.backend.cell_management.admin.views import cell_info_admin as cell_info_admin_module


# Inject our custom admin URL into the admin site's URL patterns so it lives under the
# 'admin' namespace and can be reversed by Unfold and Django admin internals.
def _inject_custom_admin_urls(orig_get_urls):
    """
    Inject custom admin URLs ahead of the default admin URLs.

    :param orig_get_urls: The original admin.site.get_urls callable

    :raise: None

    :return: A wrapper callable that returns combined URL patterns
    """

    def _get_urls():
        custom_urls = [
            path(
                "cell_management/cellinfo/",
                admin.site.admin_view(cell_info_admin_module.CellInfoAdminView.as_view()),
                name="cell_management_cellinfo_changelist",
            ),
            # Mocked filter endpoints for the Cell Info page
            path(
                "cell_management/cellinfo/filters/count/",
                admin.site.admin_view(cell_info_admin_module.cellinfo_filters_count),
                name="cell_management_cellinfo_filters_count",
            ),
            # Extract Curriculum page
            path(
                "cell_management/cellinfo/extract/",
                admin.site.admin_view(cell_info_admin_module.ExtractCurriculumAdminView.as_view()),
                name="cell_management_cellinfo_extract",
            ),
        ]
        return custom_urls + orig_get_urls()

    return _get_urls


# Only wrap once to avoid duplicating URLs on autoreload
if not getattr(admin.site, "_cellarium_custom_urls_injected", False):
    admin.site.get_urls = _inject_custom_admin_urls(admin.site.get_urls)  # type: ignore[method-assign]
    setattr(admin.site, "_cellarium_custom_urls_injected", True)

urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/", include("cellarium.nexus.backend.api_urls")),
    path("", RedirectView.as_view(url="/admin/", permanent=False)),
]
