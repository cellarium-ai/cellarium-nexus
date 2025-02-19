from django.contrib.contenttypes.models import ContentType
from django.db import IntegrityError
from django.db.models import DateTimeField
from django.utils import timezone
from nexus.backend.app_modules.cell_management import models
from rest_framework import serializers
from rest_framework.exceptions import ValidationError


class IngestInfoSerializer(serializers.ModelSerializer):
    id = serializers.UUIDField(read_only=True)
    bigquery_dataset = serializers.SlugRelatedField(slug_field="name", queryset=models.BigQueryDataset.objects.all())
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


class BulkCreateUpdateListSerializer(serializers.ListSerializer):
    def create(self, validated_data):

        result = [self.child.Meta.model(**attrs) for attrs in validated_data]

        try:
            self.child.Meta.model.objects.bulk_create(result)
        except IntegrityError as e:
            raise ValidationError(e)

        return result

    def update(self, instances, validated_data):
        instance_hash = {index: instance for index, instance in enumerate(instances)}
        result = [self.child.update(instance_hash[index], attrs) for index, attrs in enumerate(validated_data)]

        writable_fields = [x for x in self.child.Meta.fields if x not in self.child.Meta.read_only_fields]

        # Add any fields with auto_now or auto_now_add attribute set to True to writable_fields
        for field in self.child.Meta.model._meta.get_fields():
            if isinstance(field, DateTimeField) and (field.auto_now or field.auto_now_add):
                writable_fields.append(field.name)

        last_modified = timezone.now()
        for instance in result:
            # Update each readonly field with current timestamp
            for name in writable_fields:
                setattr(instance, name, last_modified)

        try:
            self.child.Meta.model.objects.bulk_update(result, writable_fields)
        except IntegrityError as e:
            raise serializers.ValidationError(e)

        return result


class CellInfoSerializer(serializers.ModelSerializer):
    id = serializers.IntegerField(required=False)
    ingest_id = serializers.PrimaryKeyRelatedField(source="ingest", queryset=models.IngestInfo.objects.all())

    class Meta:
        model = models.CellInfo
        list_serializer_class = BulkCreateUpdateListSerializer
        fields = (
            "id",
            "original_id",
            "ingest_id",
            "metadata_extra",
            "donor_id",
            "cell_type",
            "assay",
            "development_stage",
            "tissue",
            "disease",
            "organism",
            "self_reported_ethnicity",
            "sex",
            "suspension_type",
            "total_mrna_umis",
            "cell_type_ontology_term_id",
            "assay_ontology_term_id",
            "development_stage_ontology_term_id",
            "tissue_ontology_term_id",
            "disease_ontology_term_id",
            "organism_ontology_term_id",
            "self_reported_ethnicity_ontology_term_id",
            "sex_ontology_term_id",
        )


class FeatureInfoSerializer(serializers.ModelSerializer):
    id = serializers.IntegerField(required=False)
    ingest_id = serializers.PrimaryKeyRelatedField(source="ingest", queryset=models.IngestInfo.objects.all())

    class Meta:
        model = models.CellFeatureInfo
        fields = [
            "id",
            "ensemble_id",
            "symbol",
            "biotype",
            "is_filtered",
            "reference",
            "ingest_id",
            "metadata_extra",
        ]


class ReserveIndexesSerializer(serializers.Serializer):
    batch_size = serializers.IntegerField(min_value=1, write_only=True)
    index_start = serializers.IntegerField(read_only=True)
    index_end = serializers.IntegerField(read_only=True)


class IngestFromAvroSerializer(serializers.Serializer):
    stage_dir = serializers.CharField(help_text="Base staging directory path")
    ingest_nexus_uuid = serializers.UUIDField(help_text="UUID of the ingest process")

    def validate(self, data):
        """
        Validate that the ingest exists and is in the correct state.

        :param data: Serializer data

        :raise ValidationError: if validation fails

        :return: Validated data
        """
        try:
            ingest = models.IngestInfo.objects.get(nexus_uuid=data["ingest_nexus_uuid"])
            if ingest.status != models.IngestInfo.STATUS_STARTED:
                raise ValidationError(
                    f"Ingest with UUID {data['ingest_nexus_uuid']} is in {ingest.status} state. Expected {models.IngestInfo.STATUS_STARTED}"
                )
            data["ingest"] = ingest
            return data
        except models.IngestInfo.DoesNotExist:
            raise ValidationError(f"Ingest with UUID {data['ingest_nexus_uuid']} does not exist")
