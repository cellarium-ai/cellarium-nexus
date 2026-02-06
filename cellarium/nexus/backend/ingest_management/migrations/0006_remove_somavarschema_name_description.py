from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("ingest_management", "0005_somavarschema_somaingestschema_somaobscolumnschema"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="somavarschema",
            name="description",
        ),
        migrations.RemoveField(
            model_name="somavarschema",
            name="name",
        ),
    ]
