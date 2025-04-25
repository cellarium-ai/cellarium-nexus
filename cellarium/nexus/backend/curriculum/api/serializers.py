from django.contrib.auth import get_user_model
from rest_framework import serializers

from cellarium.nexus.backend.curriculum.models import Curriculum

User = get_user_model()


class CurriculumSerializer(serializers.ModelSerializer):
    """
    Serialize Curriculum model.

    :raise ValidationError: If validation fails for any field
    """

    creator_id = serializers.PrimaryKeyRelatedField(queryset=User.objects.all(), source="creator")

    class Meta:
        model = Curriculum
        fields = [
            "id",
            "name",
            "creator_id",
            "cell_count",
            "extract_bin_size",
            "extract_files_path",
            "metadata_file_path",
            "created_at",
            "filters_json",
            "status",
        ]
        read_only_fields = ["id", "created_at"]
