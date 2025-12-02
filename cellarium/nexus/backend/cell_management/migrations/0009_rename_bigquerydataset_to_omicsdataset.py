"""Rename BigQueryDataset to OmicsDataset and add backend/uri fields."""

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("cell_management", "0008_make_symbol_nullable"),
        ("ingest_management", "0003_alter_indextracking_unique_together_and_more"),
    ]

    operations = [
        # Step 1: Rename the model
        migrations.RenameModel(
            old_name="BigQueryDataset",
            new_name="OmicsDataset",
        ),
        # Step 2: Add the new 'backend' field with default 'bigquery'
        migrations.AddField(
            model_name="omicsdataset",
            name="backend",
            field=models.CharField(
                choices=[("bigquery", "BigQuery"), ("tiledb_soma", "TileDB SOMA")],
                default="bigquery",
                help_text="The storage backend for this dataset",
                max_length=32,
                verbose_name="backend",
            ),
        ),
        # Step 3: Add the new 'uri' field (nullable)
        migrations.AddField(
            model_name="omicsdataset",
            name="uri",
            field=models.CharField(
                blank=True,
                help_text="URI for the dataset (e.g., gs:// path for SOMA, or BigQuery dataset name)",
                max_length=512,
                null=True,
                unique=True,
                verbose_name="URI",
            ),
        ),
        # Step 4: Alter link field to add help_text
        migrations.AlterField(
            model_name="omicsdataset",
            name="link",
            field=models.CharField(
                blank=True,
                help_text="Link to the Google Cloud dashboard for this dataset",
                max_length=512,
                null=True,
                unique=True,
                verbose_name="link",
            ),
        ),
        # Step 5: Update verbose names
        migrations.AlterModelOptions(
            name="omicsdataset",
            options={"verbose_name": "omics dataset", "verbose_name_plural": "omics datasets"},
        ),
        # Step 6: Rename index on FeatureInfo
        migrations.RenameIndex(
            model_name="featureinfo",
            new_name="cell_manage_feature_19ee76_idx",
            old_name="cell_manage_feature_idx",
        ),
    ]
