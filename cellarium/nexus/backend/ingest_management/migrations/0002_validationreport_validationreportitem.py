# Generated by Django 5.1.4 on 2025-04-25 18:38

import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("ingest_management", "0001_initial"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name="ValidationReport",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("created_at", models.DateTimeField(auto_now_add=True, verbose_name="created at")),
                (
                    "creator",
                    models.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.SET_NULL,
                        related_name="validation_reports",
                        to=settings.AUTH_USER_MODEL,
                        verbose_name="creator",
                    ),
                ),
            ],
            options={
                "verbose_name": "validation report",
                "verbose_name_plural": "validation reports",
                "ordering": ["-created_at"],
            },
        ),
        migrations.CreateModel(
            name="ValidationReportItem",
            fields=[
                ("id", models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("input_file_gcs_path", models.CharField(max_length=1024, verbose_name="input file GCS path")),
                ("validator_name", models.CharField(max_length=255, verbose_name="validator name")),
                ("is_valid", models.BooleanField(verbose_name="is valid")),
                ("message", models.TextField(blank=True, null=True, verbose_name="message")),
                ("created_at", models.DateTimeField(auto_now_add=True, verbose_name="created at")),
                (
                    "report",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        related_name="items",
                        to="ingest_management.validationreport",
                        verbose_name="report",
                    ),
                ),
            ],
            options={
                "verbose_name": "validation report item",
                "verbose_name_plural": "validation report items",
                "ordering": ["-created_at"],
            },
        ),
    ]
