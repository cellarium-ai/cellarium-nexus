from django.urls import include, path

urlpatterns = [
    path("cell_management/", include("nexus.backend.app_modules.cell_management.api.urls")),
]
