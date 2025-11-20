from pydantic import BaseModel, Field

from cellarium.nexus.omics_datastore.bq_ops.bq_avro_schemas.custom_types import JSONBQField


class IngestInfoBQAvroSchema(BaseModel):
    id: int = Field(title="Ingest ID")
    metadata_extra: JSONBQField | None = Field(default="{}", description="Unstructured metadata, optional")


class FeatureInfoBQAvroSchema(BaseModel):
    id: int = Field(title="id", description="Primary key, unique identifier")
    ensemble_id: str = Field(title="ensemble id", description="Original identifier")
    symbol: str | None = Field(title="symbol", description="Name of the feature")
    tag: str | None = Field(default=None, title="tag")
    ingest_id: int = Field(title="ingest id", description="Foreign key referencing IngestInfo")
    metadata_extra: JSONBQField = Field(description="Extra metadata for the variable, optional")
    biotype: str | None = Field(default=None, title="biotype", description="Biotype of the feature, optional")
    is_filtered: bool | None = Field(
        default=None, title="is filtered", description="Indicates if the feature is filtered, optional"
    )
    reference: str = Field(title="reference", description="Reference information for the feature")

    class Config:
        title = "FeatureInfo"
        description = "Pydantic model for feature information, derived from Django model."
        from_attributes = True


class CellInfoBQAvroSchema(BaseModel):
    id: int = Field(title="ID")
    original_id: str = Field(title="Original ID")
    ingest_id: int = Field(title="Ingest ID")
    metadata_extra: JSONBQField = Field(title="Obs Metadata Extra")
    tag: str | None = Field(default=None, title="tag")
    # Cell Features
    donor_id: str | None = Field(default=None, title="Donor ID")
    cell_type: str | None = Field(default=None, title="Cell Type")
    assay: str | None = Field(default=None, title="Assay")
    development_stage: str | None = Field(default=None, title="Development Stage")
    tissue: str | None = Field(default=None, title="Tissue")
    disease: str | None = Field(default=None, title="Disease")
    organism: str | None = Field(default=None, title="Organism")
    self_reported_ethnicity: str | None = Field(default=None, title="Self-Reported Ethnicity")
    sex: str | None = Field(default=None, title="Sex")
    suspension_type: str | None = Field(default=None, title="Suspension Type")
    total_mrna_umis: int | None = Field(default=None, title="Total mRNA UMIs")

    # Cell Features Ontology Term IDs
    cell_type_ontology_term_id: str | None = Field(default=None, title="Cell Type Ontology Term ID")
    assay_ontology_term_id: str | None = Field(default=None, title="Assay Ontology Term ID")
    development_stage_ontology_term_id: str | None = Field(default=None, title="Development Stage Ontology Term ID")
    tissue_ontology_term_id: str | None = Field(default=None, title="Tissue Ontology Term ID")
    disease_ontology_term_id: str | None = Field(default=None, title="Disease Ontology Term ID")
    organism_ontology_term_id: str | None = Field(default=None, title="Organism Ontology Term ID")
    self_reported_ethnicity_ontology_term_id: str | None = Field(
        default=None, title="Self-Reported Ethnicity Ontology Term ID"
    )
    sex_ontology_term_id: str | None = Field(default=None, title="Sex Ontology Term ID")


class RawCountMatrixCOOBQAvroSchema(BaseModel):
    cell_id: int = Field(title="Cell ID")
    feature_id: int = Field(title="Feature ID")
    raw_count: int = Field(title="Raw Count Value")
