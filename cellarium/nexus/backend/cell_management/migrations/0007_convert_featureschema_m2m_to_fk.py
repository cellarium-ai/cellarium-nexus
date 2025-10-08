# Generated manually on 2025-10-08 13:05

import django.db.models.deletion
from django.db import migrations, models


def delete_all_features_and_schemas(apps, schema_editor):
    """
    Delete all existing FeatureInfo and FeatureSchema records.

    This is a clean slate approach - all existing data will be removed.
    """
    FeatureInfo = apps.get_model("cell_management", "FeatureInfo")
    FeatureSchema = apps.get_model("cell_management", "FeatureSchema")

    feature_count = FeatureInfo.objects.count()
    schema_count = FeatureSchema.objects.count()

    # Delete all features and schemas
    FeatureInfo.objects.all().delete()
    FeatureSchema.objects.all().delete()

    if feature_count > 0 or schema_count > 0:
        print(f"Deleted {feature_count} features and {schema_count} schemas")


def reverse_delete(apps, schema_editor):
    """
    Reverse migration: no-op since we can't restore deleted data.
    """
    print("Warning: Cannot restore deleted features and schemas")


class Migration(migrations.Migration):
    dependencies = [
        ("cell_management", "0006_remove_cellinfo_bigquery_dataset_and_more"),
    ]

    operations = [
        # Step 1: Remove M2M field from FeatureSchema (this drops the join table)
        migrations.RemoveField(
            model_name="featureschema",
            name="features",
        ),
        # Step 2: Delete all existing features and schemas (clean slate)
        migrations.RunPython(
            code=delete_all_features_and_schemas,
            reverse_code=reverse_delete,
        ),
        # Step 3: Remove old unique_together constraint on FeatureInfo
        migrations.AlterUniqueTogether(
            name="featureinfo",
            unique_together=set(),
        ),
        # Step 4: Add FK field to FeatureInfo
        migrations.AddField(
            model_name="featureinfo",
            name="feature_schema",
            field=models.ForeignKey(
                on_delete=django.db.models.deletion.CASCADE,
                related_name="features",
                to="cell_management.featureschema",
                verbose_name="feature schema",
            ),
        ),
        # Step 5: Add new unique constraint on (feature_schema, ensemble_id)
        migrations.AlterUniqueTogether(
            name="featureinfo",
            unique_together={("feature_schema", "ensemble_id")},
        ),
        # Step 6: Add index on (feature_schema, ensemble_id)
        migrations.AddIndex(
            model_name="featureinfo",
            index=models.Index(fields=["feature_schema", "ensemble_id"], name="cell_manage_feature_idx"),
        ),
    ]
