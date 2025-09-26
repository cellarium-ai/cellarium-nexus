"""
Utility functions for exporting data to CSV files.
"""

import csv
import datetime

from django.db import models
from django.http import HttpResponse


def export_model_queryset_to_csv(
    queryset: models.QuerySet,
    exclude_fields: list[str] = None,
    filename_prefix: str = None,
    identifier: str = None,
) -> HttpResponse:
    """
    Export a Django model queryset to a CSV file.

    Generate a CSV file containing all fields from the model instances in the queryset.
    The CSV file is named using the provided prefix, identifier, and current timestamp.
    If filename_prefix is not provided, the model name (lowercase) is used.
    If identifier is not provided, the current timestamp is used.

    :param queryset: Django model queryset to export
    :param exclude_fields: Optional list of field names to exclude from the export
    :param filename_prefix: Optional prefix for the CSV filename
    :param identifier: Optional identifier to include in the filename

    :raise ValueError: If queryset is empty
    :raise AttributeError: If a field cannot be accessed on a model instance

    :return: HTTP response with CSV file attachment
    """
    if not queryset.exists():
        raise ValueError("Cannot export an empty queryset")

    # Get the model class from the queryset
    model = queryset.model

    # Set default filename prefix if not provided
    if filename_prefix is None:
        filename_prefix = model._meta.model_name.lower()

    # Set default identifier if not provided
    if identifier is None:
        identifier = "export"

    # Format timestamp for filename
    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    filename = f"{filename_prefix}_{identifier}_{timestamp}.csv"

    # Create HTTP response with CSV content type
    response = HttpResponse(content_type="text/csv")
    response["Content-Disposition"] = f"attachment; filename={filename}"

    # Get all field names from the model, excluding any specified fields
    exclude_fields = exclude_fields or []
    fields = [field.name for field in model._meta.fields if field.name not in exclude_fields]

    # Create CSV writer and write header row
    writer = csv.writer(response)
    writer.writerow(fields)

    # Write data rows
    for item in queryset:
        row = []
        for field in fields:
            value = getattr(item, field)
            # Check if this is a foreign key field and get its ID instead of string representation
            field_obj = model._meta.get_field(field)
            if field_obj.is_relation and field_obj.many_to_one:
                # For foreign keys, use the ID value
                if value is not None:
                    value = value.pk
            row.append(value)
        writer.writerow(row)

    return response
