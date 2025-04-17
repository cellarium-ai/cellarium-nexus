from django.urls import include, path

urlpatterns = [
    path("ingest-management/", include("cellarium.nexus.backend.ingest_management.api.urls")),
    path("curriculum/", include("cellarium.nexus.backend.curriculum.api.urls")),
]
