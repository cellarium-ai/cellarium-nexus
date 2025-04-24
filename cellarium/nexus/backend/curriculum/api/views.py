from typing import override

from rest_framework import viewsets

from cellarium.nexus.backend.curriculum.api.serializers import CurriculumSerializer
from cellarium.nexus.backend.curriculum.models import Curriculum


class CurriculumViewSet(viewsets.ModelViewSet):
    """
    Handle CRUD operations for Curriculum model.

    List, create, retrieve, update and delete Curriculum instances.
    Creator is automatically set to the authenticated user on creation.
    """

    serializer_class = CurriculumSerializer
    queryset = Curriculum.objects.all()
    lookup_field = "name"
    lookup_url_kwarg = "name"
