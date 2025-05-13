"""
Utility functions for handling filters in the cell management app.
"""

import json
import logging
from typing import Any

from django.http import HttpRequest

from cellarium.nexus.backend.cell_management.models import BigQueryDataset

logger = logging.getLogger(__name__)


def extract_filters_from_django_admin_request(
    request: HttpRequest, dataset_filter_key: str = "bigquery_dataset__id__exact"
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
                bigquery_dataset = BigQueryDataset.objects.get(id=int(value))
                logger.info(f"Using dataset from filter: {bigquery_dataset.name}")
            except (BigQueryDataset.DoesNotExist, ValueError) as e:
                logger.error(f"Error retrieving BigQuery dataset: {str(e)}")
                raise
            continue

        # Skip parameters that don't correspond to model fields
        if key.endswith("_exclude"):
            continue

        # Get the base key and check for exclude flag
        is_exclude = False

        exclude_value = request.GET.get(f"{key}_exclude")
        if exclude_value == ["on"]:
            is_exclude = True
            logger.info(f"Found exclude flag for {key}")

        # Handle range filters
        # if key.endswith("from") or key.endswith("to"):
        if "from" in key or "to" in key:
            value = value[0]
            from_value = value if key.endswith("_from") else None
            to_value = value if key.endswith("_to") else None
            base_key = key.replace("_from", "").replace("_to", "")

            if from_value or to_value:
                if from_value:
                    try:
                        numeric_value = float(from_value)
                        filters[f"c.{base_key}__gte"] = numeric_value
                    except ValueError:
                        filters[f"c.{base_key}__gte"] = from_value

                if to_value:
                    try:
                        numeric_value = float(to_value)
                        filters[f"c.{base_key}__lte"] = numeric_value
                    except ValueError:
                        filters[f"c.{base_key}__lte"] = to_value
                continue

        # Handle multiple values
        def update_filters_with_multiple_values(
            _filters: dict[str, Any], _is_exclude: bool, _key: str, _values: list[str]
        ):
            if _values:
                _filter_key = f"c.{_key}__not_in" if _is_exclude else f"c.{_key}__in"
                _filters[_filter_key] = _values

        def update_filters_with_single_value(_filters: dict[str, Any], _is_exclude: bool, _key: str, _value: str):
            if _is_exclude:
                _filters[f"c.{_key}__not_eq"] = _value
            else:
                _filters[f"c.{_key}__eq"] = _value

        if isinstance(value, list):
            # Handle list values
            if len(value) > 1:
                values = value
                update_filters_with_multiple_values(filters, is_exclude, key, values)
            else:
                value = value[0]
                update_filters_with_single_value(filters, is_exclude, key, value)
        else:
            # Handle string value
            if "," in value:
                values = [v.strip() for v in value.split(",") if v.strip()]
                update_filters_with_multiple_values(filters, is_exclude, key, values)
            else:
                update_filters_with_single_value(filters, is_exclude, key, value)

    # Log the final filters
    logger.info(f"Extracted filters: {filters}")

    return filters, bigquery_dataset


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
