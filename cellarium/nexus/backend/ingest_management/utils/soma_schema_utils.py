import io
from dataclasses import dataclass

import pandas as pd


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
