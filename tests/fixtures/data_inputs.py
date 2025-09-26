import pathlib

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cellarium.nexus.omics_datastore.bq_ops import constants as ingest_constants


@pytest.fixture()
def small_csr_matrix() -> sp.csr_matrix:
    """
    Build a small deterministic CSR matrix for testing.

    :return: 4x3 CSR matrix with a few non-zero entries
    """
    data = np.array([1, 2, 3, 4, 5], dtype=np.int64)
    rows = np.array([0, 1, 2, 3, 3], dtype=np.int64)
    cols = np.array([0, 1, 2, 0, 2], dtype=np.int64)
    mat = sp.coo_matrix((data, (rows, cols)), shape=(4, 3)).tocsr()
    return mat


@pytest.fixture()
def small_anndata(tmp_path: pathlib.Path, small_csr_matrix: sp.csr_matrix) -> tuple[anndata.AnnData, pathlib.Path]:
    """
    Create a small AnnData with obs/var/uns and write it to an .h5ad file.

    Return both the in-memory object and the file path for read-path tests.

    :param tmp_path: Temporary directory path provided by pytest
    :param small_csr_matrix: Small deterministic CSR matrix

    :return: Tuple of (AnnData instance, h5ad file path)
    """
    n_obs, n_vars = small_csr_matrix.shape

    obs = pd.DataFrame(
        data={
            "sample": ["A", "A", "B", "B"],
            "quality": [0.1, 0.9, 0.8, 0.2],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    var = pd.DataFrame(
        data={
            "gene": ["g1", "g2", "g3"],
            "chrom": ["1", "X", "MT"],
        },
        index=[f"gene_{j}" for j in range(n_vars)],
    )

    uns = {
        "title": "test_dataset",
        "description": "unit-test",
        # large value intended to exceed small metadata limit in tests
        "large_vec": np.arange(0, 1000, dtype=np.int64),
    }

    ad = anndata.AnnData(X=small_csr_matrix, obs=obs, var=var, uns=uns)

    h5_path = tmp_path / "input.h5ad"
    ad.write_h5ad(filename=h5_path)

    return ad, h5_path


@pytest.fixture()
def obs_var_mappings() -> dict[str, dict[str, str]]:
    """
    Provide example column mappings for obs and var, including index mapping.

    :return: Mapping dictionary for obs and var
    """
    return {
        "obs_mapping": {
            "index": ingest_constants.OBS_CELL_INFO_ORIGINAL_ID,
            "sample": "sample_name",
        },
        "var_mapping": {
            "index": ingest_constants.VAR_FEATURE_INFO_ORIGINAL_ID,
            "gene": "gene_name",
        },
    }
