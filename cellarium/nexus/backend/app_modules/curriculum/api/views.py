from rest_framework import viewsets
from rest_framework.permissions import IsAuthenticated

from nexus.backend.app_modules.curriculum.api.serializers import CurriculumSerializer
from nexus.backend.app_modules.curriculum.models import Curriculum


class CurriculumViewSet(viewsets.ModelViewSet):
    """
    Handle CRUD operations for Curriculum model.

    List, create, retrieve, update and delete Curriculum instances.
    Creator is automatically set to the authenticated user on creation.
    """
    serializer_class = CurriculumSerializer
    permission_classes = [IsAuthenticated]
    queryset = Curriculum.objects.all()

    def perform_create(self, serializer: CurriculumSerializer) -> None:
        """
        Create a new Curriculum instance.

        Set the creator to the authenticated user.

        :param serializer: Validated serializer instance
        """
        serializer.save(creator=self.request.user) 