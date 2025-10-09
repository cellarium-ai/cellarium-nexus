"""
Utility functions for feature schema management.
"""

import logging
from typing import BinaryIO

import pandas as pd

from cellarium.nexus.backend.cell_management.models import FeatureInfo, FeatureSchema
from cellarium.nexus.backend.cell_management.utils.exceptions import AdminUtilsBaseException

logger = logging.getLogger(__name__)


class FeatureSchemaError(AdminUtilsBaseException):
    """
    Base exception for feature schema operations.
    """

    pass


class InvalidFeatureCSVError(FeatureSchemaError):
    """
    Raised when CSV file is invalid or missing required columns.
    """

    pass


class DuplicateSchemaNameError(FeatureSchemaError):
    """
    Raised when attempting to create a schema with an existing name.
    """

    pass


def parse_feature_csv(csv_file: BinaryIO | str) -> pd.DataFrame:
    """
    Parse and validate CSV file containing feature data.

    Validate that the CSV contains required columns: ensemble_id and symbol.
    The function intentionally does not verify that ensemble_id values are
    present in a particular dataset. Feature schemas can be reused across
    multiple dataset extracts, so we defer dataset-specific validation to the
    extract workflow (see discussion #43).

    :param csv_file: File-like object or path to CSV file

    :raise InvalidFeatureCSVError: If CSV is missing required columns or is empty
    :raise pd.errors.EmptyDataError: If CSV file is empty

    :return: Validated DataFrame with ensemble_id and symbol columns
    """
    try:
        df = pd.read_csv(csv_file)
    except pd.errors.EmptyDataError as e:
        raise InvalidFeatureCSVError("CSV file is empty") from e
    except Exception as e:
        raise InvalidFeatureCSVError(f"Failed to parse CSV file: {e}") from e

    required_columns = {"ensemble_id", "symbol"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise InvalidFeatureCSVError(
            f"CSV is missing required columns: {', '.join(sorted(missing_columns))}. "
            f"Required columns are: {', '.join(sorted(required_columns))}"
        )

    if df.empty:
        raise InvalidFeatureCSVError("CSV file contains no data rows")

    # Validate no null values in ensemble_id (symbol can be null)
    if df["ensemble_id"].isnull().any():
        raise InvalidFeatureCSVError("CSV contains null values in ensemble_id column")

    # Check for duplicate ensemble_ids
    duplicate_ensemble_ids = df[df.duplicated(subset=["ensemble_id"], keep=False)]["ensemble_id"].unique()
    if len(duplicate_ensemble_ids) > 0:
        raise InvalidFeatureCSVError(
            f"CSV contains duplicate ensemble_id values: {', '.join(duplicate_ensemble_ids.tolist())}. "
            f"Each ensemble_id must be unique within a schema."
        )

    logger.info(f"Successfully parsed CSV with {len(df)} features")
    return df


def create_feature_schema_from_dataframe(
    schema_name: str,
    features_df: pd.DataFrame,
) -> tuple[FeatureSchema, dict[str, int]]:
    """
    Create FeatureSchema and associated FeatureInfo records from DataFrame.

    Create all features as new records belonging to the schema. Features are not
    shared across schemas and will be automatically deleted when the schema is deleted.

    :param schema_name: Name for the new schema
    :param features_df: DataFrame with ensemble_id and symbol columns

    :raise DuplicateSchemaNameError: If schema with this name already exists
    :raise ValueError: If features_df is empty or missing required columns

    :return: Tuple of (created schema, stats dict with 'total_count')
    """
    # Validate schema name doesn't exist
    if FeatureSchema.objects.filter(name=schema_name).exists():
        raise DuplicateSchemaNameError(f"Schema with name '{schema_name}' already exists")

    # Validate DataFrame
    required_columns = {"ensemble_id", "symbol"}
    if not required_columns.issubset(features_df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_columns}")

    if features_df.empty:
        raise ValueError("DataFrame is empty")

    # Create the schema
    schema = FeatureSchema.objects.create(name=schema_name)
    logger.info(f"Created FeatureSchema: {schema_name}")

    # Create all features for this schema
    new_features = [
        FeatureInfo(
            ensemble_id=row["ensemble_id"],
            symbol=row["symbol"],
            feature_schema=schema,
        )
        for _, row in features_df.iterrows()
    ]

    # Bulk create all features
    FeatureInfo.objects.bulk_create(new_features)
    logger.info(f"Created {len(new_features)} FeatureInfo records for schema '{schema_name}'")

    stats = {
        "total_count": len(new_features),
    }

    return schema, stats
