# Generated by Django 5.1.4 on 2025-05-16 20:52

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("cell_management", "0004_cellinfo_bigquery_dataset"),
    ]

    operations = [
        migrations.AlterField(
            model_name="featureinfo",
            name="ensemble_id",
            field=models.CharField(max_length=255, verbose_name="ensemble id"),
        ),
        migrations.AlterUniqueTogether(
            name="featureinfo",
            unique_together={("ensemble_id", "symbol")},
        ),
    ]
