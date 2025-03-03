"""
Utility functions for handling filters in the cell management app.
"""

import json
import logging
from typing import Any

from django.http import HttpRequest
from nexus.backend.app_modules.cell_management.models import BigQueryDataset

logger = logging.getLogger(__name__)


def extract_filters_from_django_admin_request(
    request: HttpRequest, dataset_filter_key: str = "ingest__bigquery_dataset"
) -> tuple[dict[str, Any], BigQueryDataset | None]:
    """
    Extract filters from Django admin request parameters and determine the BigQuery dataset.

    Extract filters from GET parameters in the admin request and handle the BigQuery
    dataset selection based on filters or defaults.

    :param request: The HTTP request containing filter parameters
    :param dataset_filter_key: The key used for BigQuery dataset filtering

    :raise BigQueryDataset.DoesNotExist: If the selected dataset doesn't exist
    :raise ValueError: If the dataset ID is invalid

    :return: A tuple containing (filters_dict, bigquery_dataset)
    """
    filters = {}
    bigquery_dataset = None

    # Log all GET parameters for debugging
    logger.info(f"All request GET parameters: {dict(request.GET.items())}")

    # Skip processing if no GET parameters
    if not request.GET:
        logger.info("No GET parameters found in request")
        return filters, bigquery_dataset

    # Process all GET parameters to build filters
    for key, value in request.GET.items():
        # Skip empty values and pagination parameters
        if not value or key in ("p", "e", "o") or key.startswith("_"):
            logger.debug(f"Skipping parameter {key}={value}")
            continue

        # Handle the BigQuery dataset filter separately
        if key == dataset_filter_key:
            try:
                bigquery_dataset = BigQueryDataset.objects.get(id=value)
                logger.info(f"Using dataset from filter: {bigquery_dataset.name}")
            except (BigQueryDataset.DoesNotExist, ValueError) as e:
                logger.error(f"Error retrieving BigQuery dataset: {str(e)}")
                raise
            continue

        # Skip parameters that don't correspond to model fields
        if key.endswith("_exclude") or key.endswith("_from") or key.endswith("_to"):
            continue

        # Get the base key and check for exclude flag
        is_exclude = False
        base_key = key
        exclude_value = request.GET.get(f"{key}_exclude")
        if exclude_value == "on":
            is_exclude = True
            logger.info(f"Found exclude flag for {key}")

        # Handle range filters
        from_value = request.GET.get(f"{key}_from")
        to_value = request.GET.get(f"{key}_to")
        if from_value or to_value:
            if from_value:
                try:
                    numeric_value = float(from_value)
                    filters[f"c.{base_key}__gte"] = numeric_value
                    logger.info(f"Added range filter (from) for {base_key}: gte={numeric_value}")
                except ValueError:
                    filters[f"c.{base_key}__gte"] = from_value
                    logger.info(f"Added string range filter (from) for {base_key}: gte={from_value}")

            if to_value:
                try:
                    numeric_value = float(to_value)
                    filters[f"c.{base_key}__lte"] = numeric_value
                    logger.info(f"Added range filter (to) for {base_key}: lte={numeric_value}")
                except ValueError:
                    filters[f"c.{base_key}__lte"] = to_value
                    logger.info(f"Added string range filter (to) for {base_key}: lte={to_value}")
            continue

        # Handle multiple values
        if "," in value:
            values = [v.strip() for v in value.split(",") if v.strip()]
            if values:
                filter_key = f"c.{key}__not_in" if is_exclude else f"c.{key}__in"
                filters[filter_key] = values
                logger.info(f"Added multiple {'exclude' if is_exclude else 'include'} filter for {key}: {values}")
        else:
            # Handle single value
            if is_exclude:
                filters[f"c.{key}__not_eq"] = value
                logger.info(f"Added single exclude filter for {key}: {value}")
            else:
                filters[f"c.{key}__eq"] = value
                logger.info(f"Added single filter for {key}: {value}")

    # Log the final filters
    logger.info(f"Extracted filters: {filters}")

    return filters, bigquery_dataset


def get_default_dataset() -> BigQueryDataset | None:
    """
    Get the default BigQuery dataset if only one exists.

    :return: The default BigQuery dataset or None if none or multiple exist
    """
    dataset_count = BigQueryDataset.objects.count()

    if dataset_count == 1:
        bigquery_dataset = BigQueryDataset.objects.first()
        if bigquery_dataset:
            logger.info(f"Using default dataset: {bigquery_dataset.name}")
            return bigquery_dataset

    return None


def serialize_filters_to_json(filters: dict[str, Any]) -> str:
    """
    Serialize filters dictionary to a JSON string.

    :param filters: Dictionary of filters to serialize

    :return: JSON string representation of filters
    """
    return json.dumps(filters, indent=2)


def deserialize_filters_from_json(filters_json: str) -> dict[str, Any]:
    """
    Deserialize filters from a JSON string.

    :param filters_json: JSON string of filters

    :raise json.JSONDecodeError: If JSON string is invalid

    :return: Dictionary of filters
    """
    return json.loads(filters_json) if filters_json else {}
