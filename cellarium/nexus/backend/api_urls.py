from django.urls import include, path

urlpatterns = [
    path("ingest-management/", include("ingest_management.api.urls")),
    path("curriculum/", include("curriculum.api.urls")),
]
