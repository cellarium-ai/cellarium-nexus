import json
from datetime import datetime
from typing import Annotated, Any, Literal

from pydantic import BaseModel, BeforeValidator


def parse_dict(value: Any) -> dict[str, Any]:
    if isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            raise ValueError("Invalid JSON string for dictionary field")
    return value


class IngestInfoAPISchema(BaseModel):
    id: int
    nexus_uuid: str
    bigquery_dataset: str
    ingest_start_timestamp: datetime
    status: Literal["STARTED", "SUCCEEDED", "FAILED"] | None = None
    metadata_extra: Annotated[dict[str, Any] | None, BeforeValidator(parse_dict)] = None
    paper_title: str | None = None
    ingest_finish_timestamp: datetime | None = None


class FeatureInfoAPISchema(BaseModel):
    id: int
    ensemble_id: str
    symbol: str
    ingest_id: str
    biotype: str | None = None
    is_filtered: bool | None = None
    reference: str | None = None
    metadata_extra: Annotated[dict[str, Any] | None, BeforeValidator(parse_dict)] = None


class CellInfoAPISchema(BaseModel):
    id: int
    original_id: str
    ingest_id: str
    metadata_extra: Annotated[dict[str, Any] | None, BeforeValidator(parse_dict)] = None

    # Cell Features
    donor_id: str | None = None
    cell_type: str | None = None
    assay: str | None = None
    development_stage: str | None = None
    tissue: str | None = None
    disease: str | None = None
    organism: str | None = None
    self_reported_ethnicity: str | None = None
    sex: str | None = None
    suspension_type: str | None = None
    total_mrna_umis: int | None = None

    # Cell Features Ontology Term IDs
    cell_type_ontology_term_id: str | None = None
    assay_ontology_term_id: str | None = None
    development_stage_ontology_term_id: str | None = None
    tissue_ontology_term_id: str | None = None
    disease_ontology_term_id: str | None = None
    organism_ontology_term_id: str | None = None
    self_reported_ethnicity_ontology_term_id: str | None = None
    sex_ontology_term_id: str | None = None
