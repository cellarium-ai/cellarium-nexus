from django.urls import path
from nexus.backend.app_modules.cell_management.api import views

urlpatterns = [
    path("ingest/create", views.IngestCreateAPIView.as_view(), name="ingest-create"),
    path("ingest/<str:id>", views.IngestRetrieveUpdateAPIView.as_view(), name="ingest-retrieve-update"),
    path("cell-info/create-bulk", views.CellInfoListCreateAPIView.as_view(), name="cell-info-create-bulk"),
    path("cell-info/reserve-indexes", views.ReserveIndexesCellInfoAPIView.as_view(), name="cell-info-reserve-indexes"),
    path(
        "feature-info/reserve-indexes",
        views.ReserveIndexesFeatureInfoAPIView.as_view(),
        name="feature-info-reserve-indexes",
    ),
    path("feature-info/create-bulk", views.FeatureInfoListCreateAPIView.as_view(), name="feature-info-create-bulk"),
    path("ingest-from-avro/", views.IngestFromAvroView.as_view(), name="ingest-from-avro"),
]
