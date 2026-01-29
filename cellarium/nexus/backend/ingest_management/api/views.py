from rest_framework import exceptions, status
from rest_framework.generics import CreateAPIView, GenericAPIView, RetrieveUpdateAPIView
from rest_framework.request import Request
from rest_framework.response import Response
from rest_framework.views import APIView

from cellarium.nexus.backend.core.utils.reset_cache import reset_cache_and_repopulate
from cellarium.nexus.backend.ingest_management import models
from cellarium.nexus.backend.ingest_management.api import serializers
from cellarium.nexus.backend.ingest_management.services import index_tracking


class IngestCreateAPIView(CreateAPIView):
    serializer_class = serializers.IngestSerializer


class IngestRetrieveUpdateAPIView(RetrieveUpdateAPIView):
    serializer_class = serializers.IngestSerializer
    queryset = models.Ingest.objects.all()
    lookup_field = "id"
    http_method_names = ("get", "put", "patch")


class IngestFileCreateAPIView(CreateAPIView):
    serializer_class = serializers.IngestInfoSerializer


class IngestFileRetrieveUpdateAPIView(RetrieveUpdateAPIView):
    serializer_class = serializers.IngestInfoSerializer
    queryset = models.IngestInfo.objects.all()
    lookup_field = "id"
    http_method_names = ("get", "put", "patch")


class ReserveIndexesAPIViewAbstract(GenericAPIView):
    """Abstract API View to reserve indexes for a given resource key and dataset."""

    serializer_class = serializers.ReserveIndexesSerializer
    resource_key: str | None = None

    def perform_reserve(self, serializer: serializers.ReserveIndexesSerializer) -> dict[str, int]:
        """
        Handle the index reservation logic.

        :param serializer: The validated serializer instance.

        :raise exceptions.ValidationError: if resource_key is missing

        :return: Dictionary with start_index and end_index.
        """
        batch_size = serializer.validated_data["batch_size"]
        omics_dataset = serializer.validated_data["omics_dataset"]

        if not self.resource_key:
            raise exceptions.ValidationError("Resource key is not configured for this endpoint")

        start_index, end_index = index_tracking.reserve_indexes(
            omics_dataset=omics_dataset, resource_key=self.resource_key, batch_size=batch_size
        )

        return {"index_start": start_index, "index_end": end_index}

    def post(self, request: Request, *args, **kwargs) -> Response:
        """
        Handles POST requests to reserve indexes.
        """
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)

        response_data = self.perform_reserve(serializer)
        return Response(response_data)


class ReserveIndexesCellInfoAPIView(ReserveIndexesAPIViewAbstract):
    """Reserve indexes for CellInfo resource."""

    resource_key = "cell_info"


class ReserveIndexesFeatureInfoAPIView(ReserveIndexesAPIViewAbstract):
    """Reserve indexes for FeatureInfo resource."""

    resource_key = "feature_info"


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
    cached filters.
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
