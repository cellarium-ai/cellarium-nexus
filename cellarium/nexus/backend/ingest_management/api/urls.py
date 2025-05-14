from django.urls import path

from cellarium.nexus.backend.ingest_management.api import views

urlpatterns = [
    path("ingest/create", views.IngestCreateAPIView.as_view(), name="ingest-create"),
    path("ingest/<str:id>", views.IngestRetrieveUpdateAPIView.as_view(), name="ingest-retrieve-update"),
    path("ingest-from-avro/", views.IngestFromAvroView.as_view(), name="ingest-from-avro"),
    path(
        "reserve-indexes/cell-info/",
        views.ReserveIndexesCellInfoAPIView.as_view(),
        name="reserve-indexes-cell-info",
    ),
    path(
        "reserve-indexes/feature-info/",
        views.ReserveIndexesFeatureInfoAPIView.as_view(),
        name="reserve-indexes-feature-info",
    ),
    path(
        "validation-report-item/create/",
        views.ValidationReportItemCreateAPIView.as_view(),
        name="validation-report-item-create",
    ),
    path(
        "reset-cache/",
        views.ResetCacheAPIView.as_view(),
        name="reset-cache",
    ),
]
