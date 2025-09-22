from pydantic import BaseModel

from cellarium.nexus.omics_datastore.bq_avro_schemas import converter as conv_module
from cellarium.nexus.omics_datastore.bq_avro_schemas import custom_types as ct_module


class AvroModel(BaseModel):
    """
    Define a Pydantic model including Optional and JSONBQField to exercise Avro conversion.
    """

    id: int
    name: str | None
    meta: ct_module.JSONBQField
    scores: list[int]


def test_pydantic_to_avro_includes_json_field() -> None:
    """
    Convert a Pydantic model to Avro and ensure JSON fields are appended with custom sqlType.
    """
    schema = conv_module.pydantic_to_avro(pydantic_model=AvroModel)

    # basic checks about the avro schema structure
    assert schema["type"] == "record"
    assert schema["name"] == "AvroModel"
    field_names = [f["name"] for f in schema["fields"]]

    # non-JSON fields should be present
    assert "id" in field_names
    assert "name" in field_names
    assert "scores" in field_names

    # JSON field should be present with sqlType JSON
    assert "meta" in field_names
    meta_idx = field_names.index("meta")
    meta_field = schema["fields"][meta_idx]
    assert meta_field["type"]["type"] == "string"
    assert meta_field["type"]["sqlType"] == "JSON"


class BQModel(BaseModel):
    """
    Define a Pydantic model to exercise BigQuery schema conversion including lists and optionals.
    """

    id: int
    name: str | None
    meta: ct_module.JSONBQField
    tags: list[str]
    attrs: dict | None


def test_pydantic_to_bigquery_schema_fields() -> None:
    """
    Convert a Pydantic model to BigQuery schema and verify field types and modes.
    """
    fields = conv_module.pydantic_to_bigquery(pydantic_model=BQModel)

    # Map for easy lookup
    by_name = {f.name: f for f in fields}

    # id: INTEGER, REQUIRED
    assert by_name["id"].field_type == "INTEGER"
    assert by_name["id"].mode == "REQUIRED"

    # name: STRING, NULLABLE
    assert by_name["name"].field_type == "STRING"
    assert by_name["name"].mode == "NULLABLE"

    # meta: JSON, REQUIRED
    assert by_name["meta"].field_type == "JSON"
    assert by_name["meta"].mode == "REQUIRED"

    # tags: element type STRING, mode REPEATED
    assert by_name["tags"].field_type == "STRING"
    assert by_name["tags"].mode == "REPEATED"

    # attrs: JSON, NULLABLE
    assert by_name["attrs"].field_type == "JSON"
    assert by_name["attrs"].mode == "NULLABLE"
