"""
Shared fixtures for SOMA ingest integration tests.
"""

import pathlib
import shutil

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

# Apply integration marker to all tests in this directory
pytestmark = pytest.mark.integration


# =====================================================================
# Fixtures for test data
# =====================================================================


@pytest.fixture()
def sample_h5ad_paths(tmp_path: pathlib.Path) -> list[str]:
    """
    Copy sample h5ad files to temp directory and return paths.

    :param tmp_path: Temporary directory provided by pytest.

    :return: List of paths to copied h5ad files.
    """
    fixture_dir = pathlib.Path(__file__).parent.parent.parent.parent.parent / "fixtures" / "data"
    filenames = [
        "adata_ingest_test_sample-0001.h5ad",
        "adata_ingest_test_sample-0002.h5ad",
        "adata_ingest_test_sample-0003.h5ad",
    ]
    paths = []
    for fname in filenames:
        src = fixture_dir / fname
        dst = tmp_path / fname
        shutil.copyfile(src=src, dst=dst)
        paths.append(str(dst))
    return paths


@pytest.fixture()
def sample_anndata(sample_h5ad_paths: list[str]) -> anndata.AnnData:
    """
    Load the first sample h5ad file as AnnData.

    :param sample_h5ad_paths: List of h5ad file paths.

    :return: AnnData object.
    """
    return anndata.read_h5ad(filename=sample_h5ad_paths[0])


@pytest.fixture()
def var_schema_from_sample(sample_anndata: anndata.AnnData) -> ExperimentVarSchema:
    """
    Create ExperimentVarSchema from sample anndata var DataFrame.

    :param sample_anndata: Sample AnnData object.

    :return: ExperimentVarSchema with var data from sample.
    """
    return ExperimentVarSchema.from_dataframe(var_df=sample_anndata.var)


# =====================================================================
# Ingest schema fixtures
# =====================================================================


@pytest.fixture()
def ingest_schema(var_schema_from_sample: ExperimentVarSchema) -> IngestSchema:
    """
    Create an IngestSchema matching the sample h5ad test files.

    :param var_schema_from_sample: ExperimentVarSchema from sample data.

    :return: IngestSchema for validation.
    """
    return IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False),
            ObsSchemaDescriptor(name="tissue", dtype="str", nullable=True),
        ],
        var_schema=var_schema_from_sample,
        x_validation_type="count_matrix",
    )


# =====================================================================
# SOMA experiment fixtures
# =====================================================================


@pytest.fixture()
def soma_experiment_uri(tmp_path: pathlib.Path) -> str:
    """
    Provide a URI for a temporary SOMA experiment.

    :param tmp_path: Temporary directory provided by pytest.

    :return: URI string for SOMA experiment.
    """
    return str(tmp_path / "test_experiment")


# =====================================================================
# Small synthetic test data fixtures
# =====================================================================


@pytest.fixture()
def small_valid_anndata(tmp_path: pathlib.Path) -> tuple[anndata.AnnData, str]:
    """
    Create a small valid AnnData for testing validation/sanitization.

    :param tmp_path: Temporary directory provided by pytest.

    :return: Tuple of (AnnData, h5ad_path).
    """
    n_obs = 10
    n_vars = 5

    # Create sparse count matrix with valid integer values
    data = np.array([1, 2, 3, 4, 5, 1, 2, 3], dtype=np.int32)
    rows = np.array([0, 1, 2, 3, 4, 5, 6, 7])
    cols = np.array([0, 1, 2, 3, 4, 0, 1, 2])
    X = sp.csr_matrix((data, (rows, cols)), shape=(n_obs, n_vars))

    obs = pd.DataFrame(
        data={
            "cell_type": [f"type_{i % 3}" for i in range(n_obs)],
            "tissue": [f"tissue_{i % 2}" for i in range(n_obs)],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    var = pd.DataFrame(
        data={"feature_name": [f"gene_{j}" for j in range(n_vars)]},
        index=[f"ENSG{j:08d}" for j in range(n_vars)],
    )

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    h5ad_path = tmp_path / "small_valid.h5ad"
    adata.write_h5ad(filename=h5ad_path)

    return adata, str(h5ad_path)


@pytest.fixture()
def small_ingest_schema() -> IngestSchema:
    """
    Create an IngestSchema matching the small_valid_anndata fixture.

    :return: IngestSchema for small test data.
    """
    var_df = pd.DataFrame(
        data={"feature_name": [f"gene_{j}" for j in range(5)]},
        index=[f"ENSG{j:08d}" for j in range(5)],
    )
    return IngestSchema(
        obs_columns=[
            ObsSchemaDescriptor(name="cell_type", dtype="string", nullable=False),
            ObsSchemaDescriptor(name="tissue", dtype="string", nullable=True),
        ],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var_df),
        x_validation_type="count_matrix",
    )
