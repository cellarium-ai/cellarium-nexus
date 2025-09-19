from typing import Any, Optional, Type, get_args, get_origin

import py_avro_schema as pas
from google.cloud import bigquery
from pydantic import BaseModel, create_model

from cellarium.nexus.omics_datastore.bq_avro_schemas.custom_types import JSONBQField


def pydantic_to_avro(pydantic_model: Type[BaseModel]) -> dict[str, Any]:
    """
    Convert a Pydantic model class to an Avro schema.

    This function converts a Pydantic model class to an Avro schema by:
    1. Extracting fields annotated as `JSON` and handling them as special types.
    2. Using `pas.schema` to generate a base schema from the remaining fields.
    3. Adding the `JSON` fields to the Avro schema with custom attributes.

    :param pydantic_model: A class that inherits from `pydantic.BaseModel`.

    :raises ValueError: If the generated schema is invalid or unsupported.

    :return: Avro schema as a dictionary.
    """
    # Separate JSON fields and others
    json_fields: dict[str, Any] = {}
    fields_except_json: dict[str, tuple] = {}

    for field_name, field_info in pydantic_model.__annotations__.items():
        # Extract the base type and check for Optional
        origin = get_origin(field_info)
        args = get_args(field_info)

        # Determine if the field is Optional
        is_optional = origin is Optional or (len(args) == 2 and type(None) in args)
        base_type = args[0] if is_optional else field_info

        if base_type == JSONBQField:
            json_fields[field_name] = field_info
        else:
            fields_except_json[field_name] = (field_info, ...)

    # Recreate the model with the same name, excluding JSON fields
    updated_model = create_model(pydantic_model.__name__, **fields_except_json)

    try:
        # Generate Avro schema for the updated model
        avro_schema = pas.schema(py_type=updated_model)
    except Exception as e:
        raise ValueError(f"Failed to generate Avro schema: {e}")

    # Add JSON fields to the schema
    for json_field_name in json_fields:
        avro_schema["fields"].append({"name": json_field_name, "type": {"type": "string", "sqlType": "JSON"}})

    return avro_schema


PYDANTIC_TO_BIGQUERY_FIELD_MAPPING = {
    int: "INTEGER",
    float: "FLOAT",
    str: "STRING",
    bool: "BOOLEAN",
    list: "ARRAY",
    dict: "JSON",
    JSONBQField: "JSON",
}


def pydantic_to_bigquery(pydantic_model: Type[BaseModel]) -> list[bigquery.SchemaField]:
    """
    Convert a Pydantic model class to a BigQuery schema.

    This function generates a list of `bigquery.SchemaField` objects based on the annotations in the given Pydantic
    model class. The field types are mapped using a predefined type mapping (`PYDANTIC_TO_BIGQUERY_FIELD_MAPPING`), with
    "STRING" as the default type.

    :param pydantic_model: A class that inherits from :class:`pydantic.BaseModel`

    :return: A list of object :class:`bigquery.SchemaField` which describes a BigQuery table schema
    """
    schema_fields = []

    for field_name, field_info in pydantic_model.__annotations__.items():
        # Extract the base type and check for Optional
        origin = get_origin(field_info)
        args = get_args(field_info)

        # Determine if the field is Optional
        is_optional = origin is Optional or (len(args) == 2 and type(None) in args)
        base_type = args[0] if is_optional else field_info

        # Map the type to a BigQuery type
        bigquery_type = PYDANTIC_TO_BIGQUERY_FIELD_MAPPING.get(base_type, "STRING")  # Default to STRING

        # Determine the mode (REPEATED, REQUIRED, or NULLABLE)
        if origin is list:
            mode = "REPEATED"
        else:
            mode = "NULLABLE" if is_optional else "REQUIRED"

        # Create a SchemaField
        schema_fields.append(bigquery.SchemaField(name=field_name, field_type=bigquery_type, mode=mode))

    return schema_fields
