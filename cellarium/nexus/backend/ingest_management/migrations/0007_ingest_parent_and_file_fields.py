from __future__ import annotations

import django.db.models.deletion
from django.db import migrations, models


def _create_parent_ingests(apps, schema_editor) -> None:
    Ingest = apps.get_model("ingest_management", "Ingest")
    IngestInfo = apps.get_model("ingest_management", "IngestInfo")

    for ingest_info in IngestInfo.objects.all():
        ingest = Ingest.objects.create(
            omics_dataset=ingest_info.omics_dataset,
            status=ingest_info.status,
            metadata_extra=ingest_info.metadata_extra,
        )
        if ingest_info.ingest_start_timestamp:
            Ingest.objects.filter(pk=ingest.pk).update(
                ingest_start_timestamp=ingest_info.ingest_start_timestamp,
                ingest_finish_timestamp=ingest_info.ingest_finish_timestamp,
            )
        ingest_info.ingest = ingest
        ingest_info.save(update_fields=["ingest"])


class Migration(migrations.Migration):
    dependencies = [
        ("ingest_management", "0006_remove_somavarschema_name_description"),
    ]

    operations = [
        migrations.CreateModel(
            name="Ingest",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                (
                    "status",
                    models.CharField(
                        choices=[("STARTED", "Started"), ("SUCCEEDED", "Succeeded"), ("FAILED", "Failed")],
                        default="STARTED",
                        max_length=20,
                        verbose_name="status",
                    ),
                ),
                ("metadata_extra", models.JSONField(blank=True, null=True, verbose_name="metadata extra")),
                (
                    "ingest_start_timestamp",
                    models.DateTimeField(auto_now_add=True, verbose_name="ingest start timestamp"),
                ),
                (
                    "ingest_finish_timestamp",
                    models.DateTimeField(blank=True, null=True, verbose_name="ingest finish timestamp"),
                ),
                (
                    "gencode_version",
                    models.PositiveSmallIntegerField(blank=True, null=True, verbose_name="gencode version"),
                ),
                (
                    "omics_dataset",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="ingest_management_ingests",
                        to="cell_management.omicsdataset",
                        verbose_name="omics dataset",
                    ),
                ),
            ],
            options={
                "verbose_name": "ingest",
                "verbose_name_plural": "ingests",
                "ordering": ["-ingest_start_timestamp"],
            },
        ),
        migrations.AddField(
            model_name="ingestinfo",
            name="gcs_file_path",
            field=models.CharField(blank=True, max_length=1024, null=True, verbose_name="gcs file path"),
        ),
        migrations.AddField(
            model_name="ingestinfo",
            name="ingest",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="files",
                to="ingest_management.ingest",
                verbose_name="ingest",
            ),
        ),
        migrations.AddField(
            model_name="ingestinfo",
            name="tag",
            field=models.CharField(blank=True, max_length=255, null=True, verbose_name="tag"),
        ),
        migrations.AlterField(
            model_name="ingestinfo",
            name="omics_dataset",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="ingest_management_ingest_files",
                to="cell_management.omicsdataset",
                verbose_name="omics dataset",
            ),
        ),
        migrations.RunPython(_create_parent_ingests, migrations.RunPython.noop),
        migrations.AlterField(
            model_name="ingestinfo",
            name="ingest",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="files",
                to="ingest_management.ingest",
                verbose_name="ingest",
            ),
        ),
    ]
