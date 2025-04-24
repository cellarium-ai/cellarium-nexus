from django.contrib.auth import get_user_model
from rest_framework import serializers

from cellarium.nexus.backend.curriculum.models import Curriculum

User = get_user_model()


class CurriculumSerializer(serializers.ModelSerializer):
    """
    Serialize Curriculum model.

    :raise ValidationError: If validation fails for any field
    """

    creator = serializers.PrimaryKeyRelatedField(queryset=User.objects.all())

    class Meta:
        model = Curriculum
        fields = [
            "id",
            "name",
            "creator",
            "cell_count",
            "extract_bin_size",
            "extract_files_dir",
            "metadata_file_path",
            "created_at",
            "filters_json",
        ]
        read_only_fields = ["id", "created_at"]
