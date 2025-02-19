from typing import Dict

from pydantic import BaseModel, Field


class ColumnMappingSchema(BaseModel):
    """Column mapping configuration for data ingestion"""

    obs_mapping: Dict[str, str] = Field(
        default_factory=dict, description="Mapping of input obs column names to schema column names"
    )
    var_mapping: Dict[str, str] = Field(
        default_factory=dict, description="Mapping of input var column names to schema column names"
    )

    class Config:
        json_schema_extra = {
            "example": {
                "obs_mapping": {"input_cell_type": "cell_type", "input_assay": "assay"},
                "var_mapping": {
                    "input_gene_id": "ensemble_id",
                    "input_gene_name": "symbol",
                    "input_gene_type": "biotype",
                },
            }
        }
