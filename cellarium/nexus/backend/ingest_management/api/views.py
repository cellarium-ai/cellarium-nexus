from typing import Type

from django.db.models import Model

from rest_framework import status
from rest_framework.generics import CreateAPIView, GenericAPIView, RetrieveUpdateAPIView
from rest_framework.request import Request
from rest_framework.response import Response
from rest_framework.views import APIView

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.core.utils.reset_cache import reset_cache_and_repopulate
from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.api import serializers
from cellarium.nexus.backend.ingest_management.services import index_tracking


class IngestCreateAPIView(CreateAPIView):
    serializer_class = serializers.IngestInfoSerializer


class IngestRetrieveUpdateAPIView(RetrieveUpdateAPIView):
    serializer_class = serializers.IngestInfoSerializer
    queryset = models.IngestInfo.objects.all()
    lookup_field = "id"
    http_method_names = ("get", "put", "patch")


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


class ResetCacheAPIView(APIView):
    """
    API view for resetting and repopulating the cache.

    This endpoint triggers a complete cache reset and repopulation of all
    cached filters. It requires admin privileges to execute.
    """

    def post(self, request: Request) -> Response:
        """
        Handle POST request to reset and repopulate the cache.

        :param request: HTTP request

        :return: Response with the list of repopulated cache keys
        """
        repopulated_keys = reset_cache_and_repopulate()

        return Response(
            {
                "status": "success",
                "message": "Cache reset and repopulation completed successfully",
                "repopulated_keys": repopulated_keys,
                "count": len(repopulated_keys),
            },
            status=status.HTTP_200_OK,
        )
