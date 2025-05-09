from typing import Type

from django.db.models import Model
from google.api_core import exceptions as google_exceptions
from rest_framework import status
from rest_framework.generics import CreateAPIView, GenericAPIView, RetrieveUpdateAPIView
from rest_framework.request import Request
from rest_framework.response import Response
from rest_framework.views import APIView

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.api import serializers
from cellarium.nexus.backend.ingest_management.services import import_from_avro, index_tracking


class IngestCreateAPIView(CreateAPIView):
    serializer_class = serializers.IngestInfoSerializer


class IngestRetrieveUpdateAPIView(RetrieveUpdateAPIView):
    serializer_class = serializers.IngestInfoSerializer
    queryset = models.IngestInfo.objects.all()
    lookup_field = "id"
    http_method_names = ("get", "put", "patch")


class IngestFromAvroView(APIView):
    """View for ingesting CellInfo and FeatureInfo from Avro files."""

    serializer_class = serializers.IngestFromAvroSerializer

    def post(self, request):
        """
        Handle POST request to ingest Avro files.

        :param request: HTTP request

        :return: Response with ingestion results
        """
        serializer = self.serializer_class(data=request.data)
        serializer.is_valid(raise_exception=True)

        stage_dir = serializer.validated_data["stage_dir"]

        try:
            cell_info_count, feature_info_count = import_from_avro.ingest_files(stage_dir=stage_dir)

            return Response(
                {
                    "message": "Ingestion completed successfully",
                    "cell_info_count": cell_info_count,
                    "feature_info_count": feature_info_count,
                },
                status=status.HTTP_201_CREATED,
            )

        except google_exceptions.NotFound:
            return Response({"error": "One or more files not found in GCS bucket"}, status=status.HTTP_404_NOT_FOUND)
        except Exception as e:
            return Response(
                {"error": f"Error during ingestion: {str(e)}"}, status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )


class ReserveIndexesAPIViewAbstract(GenericAPIView):
    """Abstract API View to reserve indexes for a given model."""

    serializer_class = serializers.ReserveIndexesSerializer
    model_class: Type[Model] = None

    def perform_reserve(self, serializer: serializers.ReserveIndexesSerializer) -> dict[str, int]:
        """
        Handles the index reservation logic.

        :param serializer: The validated serializer instance.

        :return: Dictionary with start_index and end_index.
        """
        batch_size = serializer.validated_data["batch_size"]

        start_index, end_index = index_tracking.reserve_indexes(model_class=self.model_class, batch_size=batch_size)

        return {"start_index": start_index, "end_index": end_index}

    def post(self, request: Request, *args, **kwargs) -> Response:
        """
        Handles POST requests to reserve indexes.
        """
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)

        response_data = self.perform_reserve(serializer)

        return Response(response_data)


class ReserveIndexesCellInfoAPIView(ReserveIndexesAPIViewAbstract):
    """Reserves indexes for CellInfo model."""

    model_class = cell_models.CellInfo


class ReserveIndexesFeatureInfoAPIView(ReserveIndexesAPIViewAbstract):
    """Reserves indexes for FeatureInfo model."""

    model_class = cell_models.CellFeatureInfo


class ValidationReportItemCreateAPIView(CreateAPIView):
    """
    API view for creating ValidationReportItem entries.

    Allows creating validation report items for existing validation reports.
    """

    serializer_class = serializers.ValidationReportItemSerializer
