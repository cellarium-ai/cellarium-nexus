"""Rename bigquery_dataset to omics_dataset in IngestInfo and IndexTracking."""

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("ingest_management", "0003_alter_indextracking_unique_together_and_more"),
        ("cell_management", "0009_rename_bigquerydataset_to_omicsdataset"),
    ]

    operations = [
        # Step 1: Remove old unique_together constraint on IndexTracking
        migrations.AlterUniqueTogether(
            name="indextracking",
            unique_together=set(),
        ),
        # Step 2: Rename field in IngestInfo
        migrations.RenameField(
            model_name="ingestinfo",
            old_name="bigquery_dataset",
            new_name="omics_dataset",
        ),
        # Step 3: Rename field in IndexTracking
        migrations.RenameField(
            model_name="indextracking",
            old_name="bigquery_dataset",
            new_name="omics_dataset",
        ),
        # Step 4: Re-add unique_together constraint with new field name
        migrations.AlterUniqueTogether(
            name="indextracking",
            unique_together={("omics_dataset", "resource_key")},
        ),
        # Step 5: Alter field on IndexTracking to update verbose_name
        migrations.AlterField(
            model_name="indextracking",
            name="omics_dataset",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="ingest_management_index_trackings",
                to="cell_management.omicsdataset",
                verbose_name="omics dataset",
            ),
        ),
        # Step 6: Alter field on IngestInfo to update verbose_name
        migrations.AlterField(
            model_name="ingestinfo",
            name="omics_dataset",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="ingest_management_ingests",
                to="cell_management.omicsdataset",
                verbose_name="omics dataset",
            ),
        ),
    ]
