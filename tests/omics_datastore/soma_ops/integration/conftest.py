"""
Shared fixtures for cross ingest/extract SOMA integration tests.
"""

import pathlib

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.shared.schemas.omics_datastore import (
    ExperimentVarSchema,
    IngestSchema,
    ObsSchemaDescriptor,
)

pytestmark = pytest.mark.integration


@pytest.fixture()
def multi_ingest_anndata_inputs(tmp_path: pathlib.Path) -> dict[str, object]:
    """
    Build three AnnData objects of varying sizes with schema subsets.

    Returns a dict containing:
    - ingest_schema: IngestSchema for the full feature/obs schema
    - datasets: list of (AnnData, h5ad_path) pairs (files are not written yet)
    - expected_gene_sums: dict of feature_name -> total counts across all datasets
    - total_cells: total number of cells across all datasets
    - feature_names: ordered list of feature_name values in schema

    Note: each AnnData uses the full schema var, but only a subset of features
    have non-zero counts to simulate subset coverage.
    """
    full_features = [f"ENSG{idx:06d}" for idx in range(1, 9)]
    feature_names = [f"gene_{idx}" for idx in range(1, 9)]
    feature_name_map = dict(zip(full_features, feature_names))

    var_schema_df = pd.DataFrame({"feature_name": feature_names}, index=full_features)
    ingest_schema = IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="string", nullable=False),
            ObsSchemaDescriptor(name="tissue", dtype="string", nullable=True),
            ObsSchemaDescriptor(name="donor_id", dtype="string", nullable=True),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_schema_df),
        x_validation_type="count_matrix",
    )

    expected_gene_sums = {name: 0 for name in feature_names}
    datasets: list[tuple[anndata.AnnData, str]] = []

    def build_dataset(
        *,
        prefix: str,
        n_obs: int,
        features: list[str],
        include_tissue: bool,
        include_donor_id: bool,
        offset: int,
    ) -> None:
        n_schema_vars = len(full_features)
        n_vars = len(features)
        values = (np.arange(n_obs * n_vars, dtype=np.int32) + offset) % 7
        values = values.reshape(n_obs, n_vars)

        full_matrix = np.zeros((n_obs, n_schema_vars), dtype=np.int32)
        feature_to_idx = {feature: idx for idx, feature in enumerate(full_features)}
        for col_idx, feature in enumerate(features):
            full_matrix[:, feature_to_idx[feature]] = values[:, col_idx]

        X = sp.csr_matrix(full_matrix)

        obs_data = {"cell_type": [f"type_{i % 3}" for i in range(n_obs)]}
        if include_tissue:
            obs_data["tissue"] = [f"tissue_{i % 2}" for i in range(n_obs)]
        if include_donor_id:
            obs_data["donor_id"] = [f"donor_{i % 4}" for i in range(n_obs)]

        obs = pd.DataFrame(obs_data, index=[f"{prefix}_cell_{i}" for i in range(n_obs)])
        var = pd.DataFrame(
            {"feature_name": [feature_name_map[feature] for feature in full_features]},
            index=full_features,
        )

        adata = anndata.AnnData(X=X, obs=obs, var=var)

        # Track expected totals per feature_name
        col_sums = values.sum(axis=0)
        for feature_id, total in zip(features, col_sums):
            expected_gene_sums[feature_name_map[feature_id]] += int(total)

        h5ad_path = tmp_path / f"{prefix}.h5ad"
        datasets.append((adata, str(h5ad_path)))

    build_dataset(
        prefix="dataset_100",
        n_obs=100,
        features=full_features[:5],
        include_tissue=True,
        include_donor_id=False,
        offset=0,
    )
    build_dataset(
        prefix="dataset_150",
        n_obs=150,
        features=full_features[2:7],
        include_tissue=False,
        include_donor_id=True,
        offset=3,
    )
    build_dataset(
        prefix="dataset_200",
        n_obs=200,
        features=[full_features[i] for i in (0, 2, 4, 6)],
        include_tissue=False,
        include_donor_id=False,
        offset=5,
    )

    return {
        "ingest_schema": ingest_schema,
        "datasets": datasets,
        "expected_gene_sums": expected_gene_sums,
        "total_cells": 100 + 150 + 200,
        "feature_names": feature_names,
    }
