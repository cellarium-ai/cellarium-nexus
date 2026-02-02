import base64
import io
from dataclasses import dataclass

import pandas as pd

from cellarium.nexus.backend.ingest_management.models import IngestSchema as DjangoIngestSchema
from cellarium.nexus.shared.schemas.omics_datastore import ExperimentVarSchema, IngestSchema, ObsSchemaDescriptor


@dataclass(frozen=True)
class VarSchemaParseResult:
    dataframe: pd.DataFrame
    index_column: str


class VarSchemaCSVError(ValueError):
    pass


def parse_var_schema_csv(*, csv_file) -> VarSchemaParseResult:
    """
    Parse a var schema CSV into a DataFrame with feature IDs as index.

    The CSV must contain either a "feature_id" or "ensemble_id" column which is used as the index.
    """
    df = pd.read_csv(csv_file)
    if df.empty:
        raise VarSchemaCSVError("CSV file is empty.")

    index_column = df.columns[0]
    if df[index_column].isnull().any():
        raise VarSchemaCSVError(f"CSV index column '{index_column}' contains null values.")

    if df[index_column].duplicated().any():
        raise VarSchemaCSVError(f"CSV index column '{index_column}' contains duplicate values.")

    df = df.set_index(index_column)
    df.index = df.index.astype(str)

    return VarSchemaParseResult(dataframe=df, index_column=index_column)


def dataframe_to_parquet_bytes(*, df: pd.DataFrame) -> bytes:
    buffer = io.BytesIO()
    df.to_parquet(buffer, index=True)
    buffer.seek(0)
    return buffer.getvalue()


def django_ingest_schema_to_pydantic(*, django_schema: DjangoIngestSchema) -> IngestSchema:
    """
    Convert Django IngestSchema model to Pydantic IngestSchema.

    :param django_schema: Django IngestSchema model instance

    :return: Pydantic IngestSchema instance

    :raises AttributeError: If schema structure is invalid
    """
    obs_columns = [
        ObsSchemaDescriptor(
            name=col.name,
            dtype=col.dtype,
            nullable=col.nullable,
        )
        for col in django_schema.obs_columns.all()
    ]

    var_schema_model = django_schema.var_schema
    parquet_bytes = var_schema_model.var_parquet_file.read()
    var_parquet_b64 = base64.b64encode(parquet_bytes).decode("utf-8")
    var_schema = ExperimentVarSchema(
        _var_parquet_b64=var_parquet_b64,
        is_subset=var_schema_model.is_subset,
    )

    return IngestSchema(
        obs_columns=obs_columns,
        var_schema=var_schema,
        x_validation_type=django_schema.x_validation_type,
    )
