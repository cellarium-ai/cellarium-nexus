# Generated manually on 2025-10-08 14:36

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("cell_management", "0007_convert_featureschema_m2m_to_fk"),
    ]

    operations = [
        migrations.AlterField(
            model_name="featureinfo",
            name="symbol",
            field=models.CharField(blank=True, max_length=255, null=True, verbose_name="symbol"),
        ),
    ]
