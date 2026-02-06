# Generated migration to add measurement_name field to IngestSchema

from django.db import migrations, models
from django.utils.translation import gettext_lazy as _


class Migration(migrations.Migration):
    dependencies = [
        ("ingest_management", "0009_remove_omics_dataset_from_schemas"),
    ]

    operations = [
        migrations.AddField(
            model_name="ingestschema",
            name="measurement_name",
            field=models.CharField(
                default="RNA",
                max_length=256,
                verbose_name=_("measurement name"),
                help_text=_("Type of measurement (e.g., RNA, ATAC, OPS)"),
            ),
            preserve_default=False,
        ),
    ]
