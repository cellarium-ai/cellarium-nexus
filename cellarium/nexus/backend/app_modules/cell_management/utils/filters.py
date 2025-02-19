"""
Utility functions for handling filters in the cell management app.
"""

import logging

from django.http import HttpRequest
from nexus.backend.app_modules.cell_management.models import BigQueryDataset

logger = logging.getLogger(__name__)


def extract_filters_from_django_admin_request(
    request: HttpRequest, dataset_filter_key: str = "ingest__bigquery_dataset"
) -> tuple[dict[str, list | str], BigQueryDataset | None]:
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

    # Process all GET parameters to build filters
    for key, value in request.GET.items():
        # Skip empty values and pagination parameters
        if not value or key in ("p", "e", "o"):
            continue

        # Handle the BigQuery dataset filter separately
        if key == dataset_filter_key:
            try:
                bigquery_dataset = BigQueryDataset.objects.get(id=value)
                logger.info(f"Using dataset from filter: {bigquery_dataset.name}")
            except (BigQueryDataset.DoesNotExist, ValueError) as e:
                logger.error(f"Error retrieving BigQuery dataset: {str(e)}")
                raise
        # Add other filters to the filters dictionary
        elif not key.startswith("_"):  # Skip internal parameters
            # Handle multiple values (field__in=val1,val2,val3)
            if key.endswith("__in"):
                # Extract the field name and handle potential normalization (underscores instead of spaces)
                field_name = key.replace("__in", "")

                # Process the filter values
                processed_values = []

                if isinstance(value, str):
                    # Simple comma-separated string (most common case)
                    if "," in value:
                        # Split by comma and clean each value
                        values = [v.strip() for v in value.split(",") if v.strip()]
                        processed_values.extend(values)
                    else:
                        # Single value
                        processed_values.append(value.strip())
                elif isinstance(value, list):
                    # Already a list - Note: This might not be needed as Django typically sends strings
                    for v in value:
                        if isinstance(v, str):
                            processed_values.append(v.strip())
                        else:
                            processed_values.append(v)
                else:
                    # Single value of any type
                    processed_values.append(value)

                # Log the processed values
                logger.info(f"Processed filter values for {field_name}: {processed_values}")

                if processed_values:
                    # Format for BigQuery SQL template: field_name__in
                    filters[f"c.{field_name}__in"] = processed_values
                    logger.info(f"Added multiple filter values for {field_name}: {processed_values}")
            else:
                # Format for BigQuery SQL template: field_name__eq
                filters[f"c.{key}__eq"] = value
                logger.info(f"Added single filter value for {key}: {value}")

    # Log the filters being used
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
