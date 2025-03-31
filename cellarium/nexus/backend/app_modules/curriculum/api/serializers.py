from rest_framework import serializers

from nexus.backend.app_modules.curriculum.models import Curriculum


class CurriculumSerializer(serializers.ModelSerializer):
    """
    Serialize Curriculum model.

    :raise ValidationError: If validation fails for any field
    """

    class Meta:
        model = Curriculum
        fields = [
            "id",
            "creator",
            "cell_count",
            "extract_bin_size",
            "extract_files_dir",
            "metadata_file_path",
            "created_at",
            "filters_json",
        ]
        read_only_fields = ["id", "creator", "created_at"]
