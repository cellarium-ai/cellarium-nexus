"""
Schema definitions for feature data.
"""

from pydantic import BaseModel


class FeatureSchema(BaseModel):
    """
    Schema for feature data.

    :param id: Unique identifier for the feature
    :param symbol: Gene symbol
    :param ensemble_id: Ensemble identifier
    """

    id: int
    symbol: str
    ensemble_id: str
