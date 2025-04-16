from typing import Type

from django.db.models import Model
from google.api_core import exceptions as google_exceptions
from nexus.backend.app_modules.cell_management import models, services
from nexus.backend.app_modules.cell_management.api import serializers
from nexus.backend.app_modules.cell_management.services import import_from_avro
from rest_framework import status
from rest_framework.generics import CreateAPIView, GenericAPIView, ListCreateAPIView, RetrieveUpdateAPIView
from rest_framework.request import Request
from rest_framework.response import Response
from rest_framework.views import APIView


class IngestCreateAPIView(CreateAPIView):
    serializer_class = serializers.IngestInfoSerializer


class IngestRetrieveUpdateAPIView(RetrieveUpdateAPIView):
    serializer_class = serializers.IngestInfoSerializer
    queryset = models.IngestInfo.objects.all()
    lookup_field = "id"
    http_method_names = ("get", "put", "patch")


class CellInfoCreateAPIView(CreateAPIView):
    serializer_class = serializers.CellInfoSerializer


class BulkListCreateAPIView(ListCreateAPIView):
    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data, many=True)
        serializer.is_valid(raise_exception=True)
        self.perform_create(serializer)
        headers = self.get_success_headers(serializer.data)
        return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)


class CellInfoListCreateAPIView(BulkListCreateAPIView):
    serializer_class = serializers.CellInfoSerializer
    queryset = models.CellInfo.objects.all()


class FeatureInfoListCreateAPIView(BulkListCreateAPIView):
    serializer_class = serializers.FeatureInfoSerializer
    queryset = models.CellFeatureInfo.objects.all()


class ReserveIndexesAPIViewAbstract(GenericAPIView):
    """Abstract API View to reserve indexes for a given model."""

    serializer_class = serializers.ReserveIndexesSerializer
    cell_management_model_class: Type[Model] = None

    def perform_reserve(self, serializer: serializers.ReserveIndexesSerializer) -> dict[str, int]:
        """
        Handles the index reservation logic.

        :param serializer: The validated serializer instance.

        :return: Dictionary with start_index and end_index.
        """
        batch_size = serializer.validated_data["batch_size"]

        start_index, end_index = services.reserve_indexes(
            model_class=self.cell_management_model_class, batch_size=batch_size
        )

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

    cell_management_model_class = models.CellInfo


class ReserveIndexesFeatureInfoAPIView(ReserveIndexesAPIViewAbstract):
    """Reserves indexes for FeatureInfo model."""

    cell_management_model_class = models.CellFeatureInfo


class IngestFromAvroView(APIView):
    """View for ingesting CellInfo and FeatureInfo from Avro files."""

    serializer_class = serializers.IngestFromAvroSerializer

    def post(self, request):
        """Handle POST request to ingest Avro files.

        :param request: HTTP request
        :return: Response with ingestion results
        """
        serializer = self.serializer_class(data=request.data)
        serializer.is_valid(raise_exception=True)

        stage_dir = serializer.validated_data["stage_dir"]
        ingest = serializer.validated_data["ingest"]

        try:
            cell_info_count, feature_info_count = import_from_avro.ingest_files(stage_dir=stage_dir, ingest=ingest)

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
