import django.db.models.deletion
from django.db import migrations, models

BATCH_SIZE = 10000


def fill_bigquery_dataset(apps, schema_editor):
    CellInfo = apps.get_model("cell_management", "CellInfo")

    # Select related to avoid N+1 queries for ingest -> bigquery_dataset
    qs = CellInfo.objects.select_related("ingest").only("id", "ingest_id")

    batch = []
    for cell in qs.iterator(chunk_size=BATCH_SIZE):
        ingest = cell.ingest
        if ingest and ingest.bigquery_dataset_id:
            cell.bigquery_dataset_id = ingest.bigquery_dataset_id
            batch.append(cell)

        if len(batch) >= BATCH_SIZE:
            CellInfo.objects.bulk_update(batch, ["bigquery_dataset"])
            batch.clear()

    if batch:
        CellInfo.objects.bulk_update(batch, ["bigquery_dataset"])


class Migration(migrations.Migration):

    dependencies = [
        ("cell_management", "0003_alter_cellinfo_assay_ontology_term_id_and_more"),
    ]

    operations = [
        migrations.AddField(
            model_name="cellinfo",
            name="bigquery_dataset",
            field=models.ForeignKey(
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="cells",
                to="cell_management.bigquerydataset",
                verbose_name="BigQuery dataset",
            ),
        ),
        migrations.RunPython(fill_bigquery_dataset),
        migrations.AlterField(
            model_name="cellinfo",
            name="bigquery_dataset",
            field=models.ForeignKey(
                to="cell_management.bigquerydataset",
                null=False,
                on_delete=django.db.models.deletion.CASCADE,
                editable=False,
                verbose_name="bigquery dataset",
            ),
        ),
        migrations.AddIndex(
            model_name="cellinfo",
            index=models.Index(fields=["bigquery_dataset"], name="cellinfo_bqds_idx"),
        ),
    ]
