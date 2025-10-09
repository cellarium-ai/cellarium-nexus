from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.utils import feature_schema_utils


@pytest.fixture()
def sample_dataframe() -> pd.DataFrame:
    """
    Create a valid DataFrame for testing.

    :return: DataFrame with ensemble_id and symbol columns
    """
    return pd.DataFrame(
        {
            "ensemble_id": ["ENSG00000141510", "ENSG00000157764", "ENSG00000146648"],
            "symbol": ["TP53", "BRAF", "EGFR"],
        }
    )


# === parse_feature_csv() tests ===


def test_parse_feature_csv_valid_file(sample_csv_file: Path) -> None:
    """
    Parse valid CSV returns DataFrame with correct columns and data.

    :param sample_csv_file: Path to valid CSV file
    """
    df = feature_schema_utils.parse_feature_csv(csv_file=sample_csv_file)

    assert "ensemble_id" in df.columns
    assert "symbol" in df.columns
    assert len(df) == 3
    assert df.iloc[0]["ensemble_id"] == "ENSG00000141510"
    assert df.iloc[0]["symbol"] == "TP53"
    assert df.iloc[1]["ensemble_id"] == "ENSG00000157764"
    assert df.iloc[1]["symbol"] == "BRAF"


def test_parse_feature_csv_missing_ensemble_id_column(tmp_path: Path) -> None:
    """
    Parse CSV missing ensemble_id column raises InvalidFeatureCSVError.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "invalid.csv"
    csv_path.write_text("symbol\nTP53\nBRAF\n")

    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError) as exc_info:
        feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    assert "ensemble_id" in str(exc_info.value)
    assert "required columns" in str(exc_info.value).lower()


def test_parse_feature_csv_missing_symbol_column(tmp_path: Path) -> None:
    """
    Parse CSV missing symbol column raises InvalidFeatureCSVError.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "invalid.csv"
    csv_path.write_text("ensemble_id\nENSG00000141510\n")

    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError) as exc_info:
        feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    assert "symbol" in str(exc_info.value)
    assert "required columns" in str(exc_info.value).lower()


def test_parse_feature_csv_missing_both_columns(invalid_csv_file: Path) -> None:
    """
    Parse CSV missing both required columns raises InvalidFeatureCSVError.

    :param invalid_csv_file: Path to invalid CSV file
    """
    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError) as exc_info:
        feature_schema_utils.parse_feature_csv(csv_file=invalid_csv_file)

    error_msg = str(exc_info.value).lower()
    assert "ensemble_id" in error_msg
    assert "symbol" in error_msg
    assert "required columns" in error_msg


def test_parse_feature_csv_empty_file(tmp_path: Path) -> None:
    """
    Parse empty CSV (no data rows) raises InvalidFeatureCSVError.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "empty.csv"
    csv_path.write_text("ensemble_id,symbol\n")

    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError) as exc_info:
        feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    assert "no data rows" in str(exc_info.value).lower()


def test_parse_feature_csv_null_ensemble_id(tmp_path: Path) -> None:
    """
    Parse CSV with null ensemble_id raises InvalidFeatureCSVError.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "null_ensemble.csv"
    csv_path.write_text("ensemble_id,symbol\n,TP53\nENSG00000157764,BRAF\n")

    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError) as exc_info:
        feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    assert "null values" in str(exc_info.value).lower()


def test_parse_feature_csv_null_symbol_allowed(tmp_path: Path) -> None:
    """
    Parse CSV with null symbol is allowed (symbol is optional).

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "null_symbol.csv"
    csv_path.write_text("ensemble_id,symbol\nENSG00000141510,\nENSG00000157764,BRAF\n")

    df = feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    assert len(df) == 2
    assert pd.isna(df.iloc[0]["symbol"])
    assert df.iloc[1]["symbol"] == "BRAF"


def test_parse_feature_csv_extra_columns_preserved(tmp_path: Path) -> None:
    """
    Parse CSV with extra columns succeeds and preserves all columns.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "extra_columns.csv"
    csv_path.write_text(
        "ensemble_id,symbol,extra_col1,extra_col2\n"
        "ENSG00000141510,TP53,value1,value2\n"
        "ENSG00000157764,BRAF,value3,value4\n"
    )

    df = feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    assert "ensemble_id" in df.columns
    assert "symbol" in df.columns
    assert "extra_col1" in df.columns
    assert "extra_col2" in df.columns
    assert len(df) == 2


def test_parse_feature_csv_malformed_file(tmp_path: Path) -> None:
    """
    Parse malformed CSV raises InvalidFeatureCSVError.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "malformed.csv"
    csv_path.write_text("not a valid csv format\nrandom text\n")

    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError):
        feature_schema_utils.parse_feature_csv(csv_file=csv_path)


# === create_feature_schema_from_dataframe() tests ===


def test_create_feature_schema_from_dataframe_success(sample_dataframe: pd.DataFrame) -> None:
    """
    Create schema with features from DataFrame succeeds.

    :param sample_dataframe: Valid DataFrame with features
    """
    schema, stats = feature_schema_utils.create_feature_schema_from_dataframe(
        schema_name="test-schema",
        features_df=sample_dataframe,
    )

    assert schema.name == "test-schema"
    assert schema.features.count() == 3
    assert stats["total_count"] == 3

    # Verify features are created correctly
    feature = schema.features.get(ensemble_id="ENSG00000141510")
    assert feature.symbol == "TP53"
    assert feature.feature_schema == schema

    # Verify all features
    features = list(schema.features.all().order_by("symbol"))
    assert features[0].symbol == "BRAF"
    assert features[1].symbol == "EGFR"
    assert features[2].symbol == "TP53"


def test_create_feature_schema_from_dataframe_duplicate_schema_name(
    sample_dataframe: pd.DataFrame,
) -> None:
    """
    Create schema with existing name raises DuplicateSchemaNameError.

    :param sample_dataframe: Valid DataFrame with features
    """
    # Create first schema
    models.FeatureSchema.objects.create(name="duplicate-name")

    # Attempt to create second with same name
    with pytest.raises(feature_schema_utils.DuplicateSchemaNameError) as exc_info:
        feature_schema_utils.create_feature_schema_from_dataframe(
            schema_name="duplicate-name",
            features_df=sample_dataframe,
        )

    assert "duplicate-name" in str(exc_info.value)
    assert "already exists" in str(exc_info.value).lower()


def test_create_feature_schema_from_dataframe_empty_dataframe() -> None:
    """
    Create schema from empty DataFrame raises ValueError.
    """
    empty_df = pd.DataFrame({"ensemble_id": [], "symbol": []})

    with pytest.raises(ValueError) as exc_info:
        feature_schema_utils.create_feature_schema_from_dataframe(
            schema_name="test-schema",
            features_df=empty_df,
        )

    assert "empty" in str(exc_info.value).lower()


def test_create_feature_schema_from_dataframe_missing_ensemble_id_column() -> None:
    """
    Create schema from DataFrame missing ensemble_id raises ValueError.
    """
    invalid_df = pd.DataFrame({"symbol": ["TP53", "BRAF"]})

    with pytest.raises(ValueError) as exc_info:
        feature_schema_utils.create_feature_schema_from_dataframe(
            schema_name="test-schema",
            features_df=invalid_df,
        )

    assert "ensemble_id" in str(exc_info.value)


def test_create_feature_schema_from_dataframe_missing_symbol_column() -> None:
    """
    Create schema from DataFrame missing symbol raises ValueError.
    """
    invalid_df = pd.DataFrame({"ensemble_id": ["ENSG00000141510", "ENSG00000157764"]})

    with pytest.raises(ValueError) as exc_info:
        feature_schema_utils.create_feature_schema_from_dataframe(
            schema_name="test-schema",
            features_df=invalid_df,
        )

    assert "symbol" in str(exc_info.value)


def test_parse_feature_csv_duplicate_ensemble_ids(tmp_path: Path) -> None:
    """
    Parse CSV with duplicate ensemble_id values raises InvalidFeatureCSVError.

    :param tmp_path: Temporary directory path
    """
    csv_path = tmp_path / "duplicates.csv"
    csv_path.write_text(
        "ensemble_id,symbol\n"
        "ENSG00000141510,TP53\n"
        "ENSG00000141510,TP53_ALT\n"  # Same ensemble_id, different symbol
        "ENSG00000157764,BRAF\n"
    )

    with pytest.raises(feature_schema_utils.InvalidFeatureCSVError) as exc_info:
        feature_schema_utils.parse_feature_csv(csv_file=csv_path)

    error_msg = str(exc_info.value)
    assert "duplicate ensemble_id" in error_msg.lower()
    assert "ENSG00000141510" in error_msg
    assert "unique within a schema" in error_msg.lower()


def test_create_feature_schema_from_dataframe_cascade_delete(
    sample_dataframe: pd.DataFrame,
) -> None:
    """
    Deleting schema cascades to delete all features.

    :param sample_dataframe: Valid DataFrame with features
    """
    schema, _ = feature_schema_utils.create_feature_schema_from_dataframe(
        schema_name="test-schema",
        features_df=sample_dataframe,
    )

    feature_ids = list(schema.features.values_list("pk", flat=True))
    assert len(feature_ids) == 3

    # Delete schema
    schema.delete()

    # Verify features are deleted (CASCADE)
    assert models.FeatureInfo.objects.filter(pk__in=feature_ids).count() == 0


def test_create_feature_schema_from_dataframe_features_isolated_per_schema() -> None:
    """
    Features from different schemas are isolated (not shared).
    """
    df1 = pd.DataFrame(
        {
            "ensemble_id": ["ENSG00000141510"],
            "symbol": ["TP53"],
        }
    )
    df2 = pd.DataFrame(
        {
            "ensemble_id": ["ENSG00000141510"],
            "symbol": ["TP53"],
        }
    )

    schema1, _ = feature_schema_utils.create_feature_schema_from_dataframe(
        schema_name="schema-1",
        features_df=df1,
    )
    schema2, _ = feature_schema_utils.create_feature_schema_from_dataframe(
        schema_name="schema-2",
        features_df=df2,
    )

    # Each schema should have its own TP53 feature
    assert schema1.features.count() == 1
    assert schema2.features.count() == 1

    # They should be different objects
    feature1 = schema1.features.first()
    feature2 = schema2.features.first()
    assert feature1.pk != feature2.pk
    assert feature1.feature_schema == schema1
    assert feature2.feature_schema == schema2

    # Deleting schema1 should not affect schema2's features
    schema1.delete()
    assert models.FeatureInfo.objects.filter(pk=feature2.pk).exists()
    assert schema2.features.count() == 1


def test_create_feature_schema_from_dataframe_large_dataset() -> None:
    """
    Create schema with large DataFrame uses bulk operations efficiently.
    """
    # Create a large DataFrame
    large_df = pd.DataFrame(
        {
            "ensemble_id": [f"ENSG{i:011d}" for i in range(1000)],
            "symbol": [f"GENE{i}" for i in range(1000)],
        }
    )

    schema, stats = feature_schema_utils.create_feature_schema_from_dataframe(
        schema_name="large-schema",
        features_df=large_df,
    )

    assert schema.features.count() == 1000
    assert stats["total_count"] == 1000

    # Verify a few random features
    assert schema.features.filter(ensemble_id="ENSG00000000000").exists()
    assert schema.features.filter(symbol="GENE500").exists()
