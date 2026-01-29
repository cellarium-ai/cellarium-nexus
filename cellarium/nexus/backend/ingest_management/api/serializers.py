from django.utils import timezone
from rest_framework import serializers
from rest_framework.exceptions import ValidationError

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.ingest_management import models


def _refresh_ingest_status(ingest: models.Ingest) -> None:
    statuses = set(ingest.files.values_list("status", flat=True))
    if not statuses:
        return

    if models.IngestInfo.STATUS_FAILED in statuses:
        new_status = models.Ingest.STATUS_FAILED
    elif statuses == {models.IngestInfo.STATUS_SUCCEEDED}:
        new_status = models.Ingest.STATUS_SUCCEEDED
    else:
        new_status = models.Ingest.STATUS_STARTED

    updates: dict[str, object] = {}
    if ingest.status != new_status:
        updates["status"] = new_status

    if new_status in (models.Ingest.STATUS_FAILED, models.Ingest.STATUS_SUCCEEDED):
        if ingest.ingest_finish_timestamp is None:
            updates["ingest_finish_timestamp"] = timezone.now()
    elif ingest.ingest_finish_timestamp is not None:
        updates["ingest_finish_timestamp"] = None

    if updates:
        models.Ingest.objects.filter(pk=ingest.pk).update(**updates)


class IngestSerializer(serializers.ModelSerializer):
    id = serializers.IntegerField(read_only=True)
    omics_dataset = serializers.SlugRelatedField(
        slug_field="name",
        queryset=cell_models.OmicsDataset.objects.all(),
        required=False,
        allow_null=True,
    )
    ingest_start_timestamp = serializers.DateTimeField(read_only=True)
    status = serializers.ChoiceField(choices=models.Ingest.STATUS_CHOICES, required=False)

    class Meta:
        model = models.Ingest
        fields = (
            "id",
            "omics_dataset",
            "status",
            "metadata_extra",
            "gencode_version",
            "ingest_start_timestamp",
            "ingest_finish_timestamp",
        )


class IngestInfoSerializer(serializers.ModelSerializer):
    id = serializers.IntegerField(read_only=True)
    ingest_id = serializers.PrimaryKeyRelatedField(
        source="ingest",
        queryset=models.Ingest.objects.all(),
        required=False,
        allow_null=True,
    )
    omics_dataset = serializers.SlugRelatedField(
        slug_field="name",
        queryset=cell_models.OmicsDataset.objects.all(),
        required=False,
        allow_null=True,
    )
    ingest_start_timestamp = serializers.DateTimeField(read_only=True)
    status = serializers.ChoiceField(choices=models.IngestInfo.STATUS_CHOICES, required=False)

    class Meta:
        model = models.IngestInfo
        fields = (
            "id",
            "ingest_id",
            "nexus_uuid",
            "omics_dataset",
            "status",
            "gcs_file_path",
            "tag",
            "metadata_extra",
            "ingest_start_timestamp",
            "ingest_finish_timestamp",
        )

    def create(self, validated_data):
        ingest = validated_data.get("ingest")
        omics_dataset = validated_data.get("omics_dataset")

        if ingest is None:
            if omics_dataset is None:
                raise ValidationError("omics_dataset is required when ingest_id is not provided.")
            ingest = models.Ingest.objects.create(
                omics_dataset=omics_dataset,
                status=models.Ingest.STATUS_STARTED,
            )
            validated_data["ingest"] = ingest
        else:
            validated_data["omics_dataset"] = ingest.omics_dataset

        return super().create(validated_data)

    def update(self, instance, validated_data):
        ingest = instance.ingest
        validated_data.pop("ingest", None)
        updated = super().update(instance, validated_data)
        _refresh_ingest_status(ingest)
        return updated


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
    """Serialize input and output for reserving indexes."""

    omics_dataset = serializers.SlugRelatedField(
        slug_field="name",
        queryset=cell_models.OmicsDataset.objects.all(),
        write_only=True,
        help_text="Name of the omics dataset to scope the reservation",
    )
    batch_size = serializers.IntegerField(min_value=1, write_only=True)
    index_start = serializers.IntegerField(read_only=True)
    index_end = serializers.IntegerField(read_only=True)


class ValidationReportItemSerializer(serializers.ModelSerializer):
    """
    Serializer for ValidationReportItem model.

    Used for creating and retrieving validation report items.
    """

    report_id = serializers.PrimaryKeyRelatedField(
        source="report",
        queryset=models.ValidationReport.objects.all(),
        required=True,
        help_text="ID of the existing ValidationReport",
    )

    class Meta:
        model = models.ValidationReportItem
        fields = (
            "id",
            "report_id",
            "input_file_gcs_path",
            "validator_name",
            "is_valid",
            "message",
            "created_at",
        )
        read_only_fields = ("id", "created_at")
