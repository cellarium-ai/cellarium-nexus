from rest_framework import serializers
from rest_framework.exceptions import ValidationError

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.ingest_management import models


class IngestInfoSerializer(serializers.ModelSerializer):
    id = serializers.UUIDField(read_only=True)
    bigquery_dataset = serializers.SlugRelatedField(
        slug_field="name", queryset=cell_models.BigQueryDataset.objects.all()
    )
    ingest_start_timestamp = serializers.DateTimeField(read_only=True)
    status = serializers.ChoiceField(choices=models.IngestInfo.STATUS_CHOICES, required=False)

    class Meta:
        model = models.IngestInfo
        fields = (
            "id",
            "nexus_uuid",
            "bigquery_dataset",
            "status",
            "metadata_extra",
            "ingest_start_timestamp",
            "ingest_finish_timestamp",
        )


class IngestFromAvroSerializer(serializers.Serializer):
    stage_dir = serializers.CharField(help_text="Base staging directory path")
    ingest_id = serializers.IntegerField(help_text="ID of the ingest process")

    def validate(self, data):
        """
        Validate that the ingest exists and is in the correct state.

        :param data: Serializer data

        :raise ValidationError: if validation fails

        :return: Validated data
        """
        try:
            ingest = models.IngestInfo.objects.get(id=data["ingest_id"])
            if ingest.status != models.IngestInfo.STATUS_STARTED:
                raise ValidationError(
                    f"Ingest with ID {data['ingest_id']} is in {ingest.status} state. Expected {models.IngestInfo.STATUS_STARTED}"
                )
            data["ingest"] = ingest
            return data
        except models.IngestInfo.DoesNotExist:
            raise ValidationError(f"Ingest with ID {data['ingest_id']} does not exist")


class ReserveIndexesSerializer(serializers.Serializer):
    """Serializer for reserving indexes."""

    batch_size = serializers.IntegerField(min_value=1, write_only=True)
    index_start = serializers.IntegerField(read_only=True)
    index_end = serializers.IntegerField(read_only=True)
