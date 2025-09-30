import json

import numpy as np
import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.bq_avro_schemas.cell_management import (
    CellInfoBQAvroSchema,
    FeatureInfoBQAvroSchema,
)
from cellarium.nexus.omics_datastore.bq_ops import exceptions as ingest_exceptions
from cellarium.nexus.omics_datastore.bq_ops.ingest import create_ingest_files


def test_apply_column_mapping_regular_and_index() -> None:
    """
    Validate that _apply_column_mapping renames columns and maps index to a new column.
    """
    df = pd.DataFrame(
        data={"a": [1, 2], "b": [3, 4]},
        index=["row0", "row1"],
    )
    mapping = {"index": "row_id", "a": "alpha"}

    out = create_ingest_files._apply_column_mapping(df=df, mapping=mapping)

    assert "alpha" in out.columns
    assert "b" in out.columns
    assert "row_id" in out.columns
    assert list(out["alpha"]) == [1, 2]
    assert list(out["row_id"]) == ["row0", "row1"]


def test_ensure_schema_fields_adds_missing() -> None:
    """
    Validate that _ensure_schema_fields adds missing columns with NA values.
    """
    df = pd.DataFrame(data={"x": [1, 2]})
    schema_fields = ["x", "y", "z"]

    out = create_ingest_files._ensure_schema_fields(df=df, schema_field_names=schema_fields)

    assert set(schema_fields).issubset(set(out.columns))
    assert out["y"].isna().all()
    assert out["z"].isna().all()


def test_numpy_json_encoder_encodes_numpy_types() -> None:
    """
    Validate that NumpyJSONEncoder converts numpy types to JSON-friendly forms.
    """
    payload = {
        "i": np.int64(5),
        "f": np.float64(3.14),
        "a": np.array([1, 2, 3], dtype=np.int64),
    }

    encoded = json.dumps(obj=payload, cls=create_ingest_files.NumpyJSONEncoder)
    decoded = json.loads(encoded)

    assert decoded == {"i": 5, "f": 3.14, "a": [1, 2, 3]}


class _Unserializable:
    pass


def test_numpy_json_encoder_handles_unserializable() -> None:
    """
    Validate that NumpyJSONEncoder emits type description for unserializable objects.
    """
    payload = {"obj": _Unserializable()}
    encoded = json.dumps(obj=payload, cls=create_ingest_files.NumpyJSONEncoder)
    decoded = json.loads(encoded)
    assert "Error while encoding type _Unserializable while ingest" in decoded["obj"]


def test_required_non_optional_string_fields_feature_info() -> None:
    """
    Validate detection of required, non-optional string fields for FeatureInfo schema.
    """
    fields = create_ingest_files._required_non_optional_string_fields(model=FeatureInfoBQAvroSchema)
    # `ensemble_id` and `reference` are required strings; `symbol` is Optional[str]
    assert fields == {"ensemble_id", "reference"}


def test_required_non_optional_string_fields_cell_info() -> None:
    """
    Validate detection of required, non-optional string fields for CellInfo schema.
    """
    fields = create_ingest_files._required_non_optional_string_fields(model=CellInfoBQAvroSchema)
    # `original_id` is required str; most others are Optional[str] or non-str
    assert fields == {"original_id"}


def test_validate_required_string_columns_happy_path() -> None:
    """
    Validate that validator passes when required string columns exist and are strings.
    """
    df = pd.DataFrame(
        data={
            "ensemble_id": ["g1", "g2"],
            "reference": ["ref1", "ref2"],
        }
    )
    required = {"ensemble_id", "reference"}
    # Should not raise
    create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")


def test_validate_required_string_columns_missing_column_raises() -> None:
    """
    Validate that validator raises when required string columns are missing.
    """
    df = pd.DataFrame(data={"ensemble_id": ["g1", "g2"]})
    required = {"ensemble_id", "reference"}
    with pytest.raises(ingest_exceptions.DataValidationError) as exc:
        create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")
    assert "Missing required var string columns" in str(exc.value)


def test_validate_required_string_columns_null_values_raise() -> None:
    """
    Validate that validator raises when required string columns contain nulls.
    """
    df = pd.DataFrame(data={"ensemble_id": ["g1", None], "reference": ["ref1", "ref2"]})
    required = {"ensemble_id", "reference"}
    with pytest.raises(ingest_exceptions.DataValidationError) as exc:
        create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")
    assert "columns with null values" in str(exc.value)


def test_validate_required_string_columns_non_string_values_raise() -> None:
    """
    Validate that validator raises when required string columns contain non-string values.
    """
    df = pd.DataFrame(data={"ensemble_id": ["g1", 2], "reference": ["ref1", "ref2"]})
    required = {"ensemble_id", "reference"}
    with pytest.raises(ingest_exceptions.DataValidationError) as exc:
        create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")
    assert "columns with non-string values" in str(exc.value)


def test_validate_required_string_columns_categorical_strings_pass() -> None:
    """
    Validate that validator works when columns are pandas Categorical of strings.
    """
    df = pd.DataFrame(
        data={
            "ensemble_id": pd.Categorical(["g1", "g2"]),
            "reference": pd.Categorical(["ref1", "ref2"]),
        }
    )
    required = {"ensemble_id", "reference"}
    # Should not raise
    create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")


def test_validate_required_string_columns_categorical_non_string_values_raise() -> None:
    """
    Validate that validator raises when Categorical contains non-strings.
    """
    df = pd.DataFrame(
        data={
            "ensemble_id": pd.Categorical(["g1", 2]),
            "reference": pd.Categorical(["ref1", "ref2"]),
        }
    )
    required = {"ensemble_id", "reference"}
    with pytest.raises(ingest_exceptions.DataValidationError) as exc:
        create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")
    assert "columns with non-string values" in str(exc.value)


def test_validate_required_string_columns_pandas_stringdtype_pass() -> None:
    """
    Validate that validator works when columns use pandas StringDtype.
    """
    df = pd.DataFrame(
        data={
            "ensemble_id": pd.Series(["g1", "g2"], dtype=pd.StringDtype()),
            "reference": pd.Series(["ref1", "ref2"], dtype=pd.StringDtype()),
        }
    )
    required = {"ensemble_id", "reference"}
    # Should not raise
    create_ingest_files._validate_required_string_columns(df=df, required_cols=required, context="var")
