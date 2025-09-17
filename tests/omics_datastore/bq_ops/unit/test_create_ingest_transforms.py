import numpy as np
import pandas as pd

from cellarium.nexus.omics_datastore.bq_ops import constants as ingest_constants
from cellarium.nexus.omics_datastore.bq_ops.ingest import create_ingest_files


def _expected_total_mrna_umis(x) -> list[int]:
    sums = x.sum(axis=1)
    # handle sparse matrix .A1 vs ndarray
    if hasattr(sums, "A1"):
        return list(sums.A1.astype(int))
    return list(np.asarray(sums).astype(int))


def test_process_cell_info_obs_with_and_without_mapping(small_csr_matrix) -> None:
    """
    Validate that _process_cell_info_obs adds id/ingest_id/tag and computes total_mrna_umis.

    Also validate mapping behavior and schema/metadata partitioning.

    :raise: None
    :return: None
    """
    n_obs, _ = small_csr_matrix.shape
    obs = pd.DataFrame(
        data={
            "sample": ["A", "A", "B", "B"],
            "quality": [0.1, 0.9, 0.8, 0.2],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    # var is not used directly here
    adata = create_ingest_files.AnnData(X=small_csr_matrix, obs=obs)

    # without mapping
    df_schema, df_meta = create_ingest_files._process_cell_info_obs(
        adata=adata,
        tag="t1",
        ingest_id=123,
        start_index=1,
        end_index=n_obs,
        column_mapping=None,
    )

    assert ingest_constants.OBS_NEXUS_ID in df_schema.columns
    assert ingest_constants.OBS_INGEST_ID in df_schema.columns
    assert ingest_constants.OBS_TAG in df_schema.columns
    assert ingest_constants.OBS_TOTAL_MRNA_UMIS in df_schema.columns
    assert df_schema[ingest_constants.OBS_TOTAL_MRNA_UMIS].tolist() == _expected_total_mrna_umis(adata.X)

    # with mapping including index
    obs_map = {"index": ingest_constants.OBS_CELL_INFO_ORIGINAL_ID, "sample": "sample_name"}
    df_schema2, df_meta2 = create_ingest_files._process_cell_info_obs(
        adata=adata,
        tag="t2",
        ingest_id=456,
        start_index=10,
        end_index=9 + n_obs,
        column_mapping=obs_map,
    )

    assert ingest_constants.OBS_CELL_INFO_ORIGINAL_ID in pd.concat([df_schema2, df_meta2], axis=1).columns
    assert "sample_name" in pd.concat([df_schema2, df_meta2], axis=1).columns
    assert df_schema2[ingest_constants.OBS_TOTAL_MRNA_UMIS].tolist() == _expected_total_mrna_umis(adata.X)


def test_process_feature_info_var_with_and_without_mapping() -> None:
    """
    Validate that _process_feature_info_var adds id/ingest_id/tag and partitions columns.

    :raise: None
    :return: None
    """
    df_var = pd.DataFrame(
        data={
            "gene": ["g1", "g2", "g3"],
            "chrom": ["1", "X", "MT"],
        }
    )

    df_schema, df_meta = create_ingest_files._process_feature_info_var(
        df=df_var.copy(),
        tag="t1",
        ingest_id=99,
        start_index=1,
        end_index=3,
        column_mapping=None,
    )

    assert ingest_constants.VAR_NEXUS_ID in df_schema.columns
    assert ingest_constants.VAR_INGEST_ID in df_schema.columns
    assert ingest_constants.VAR_TAG in df_schema.columns

    var_map = {"index": ingest_constants.VAR_FEATURE_INFO_ORIGINAL_ID, "gene": "gene_name"}

    df_schema2, df_meta2 = create_ingest_files._process_feature_info_var(
        df=df_var.copy(),
        tag="t2",
        ingest_id=77,
        start_index=10,
        end_index=12,
        column_mapping=var_map,
    )

    full2 = pd.concat([df_schema2, df_meta2], axis=1)
    assert ingest_constants.VAR_FEATURE_INFO_ORIGINAL_ID in full2.columns
    assert "gene_name" in full2.columns
