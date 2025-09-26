import numpy as np
import pandas as pd
import pytest

from cellarium.nexus.omics_datastore.bq_ops import constants as ingest_constants
from cellarium.nexus.omics_datastore.bq_ops.ingest import create_ingest_files


def _expected_total_mrna_umis(x) -> list[int]:
    sums = x.sum(axis=1)
    # handle sparse matrix .A1 vs ndarray
    if hasattr(sums, "A1"):
        return list(sums.A1.astype(int))
    return list(np.asarray(sums).astype(int))


@pytest.mark.parametrize(
    "column_mapping, tag, ingest_id, start_index, expected_extra_columns",
    [
        pytest.param(None, "t1", 123, 1, [], id="without-mapping"),
        pytest.param(
            {"index": ingest_constants.OBS_CELL_INFO_ORIGINAL_ID, "sample": "sample_name"},
            "t2",
            456,
            10,
            [ingest_constants.OBS_CELL_INFO_ORIGINAL_ID, "sample_name"],
            id="with-mapping",
        ),
    ],
)
def test_process_cell_info_obs_mapping_variants(
    small_csr_matrix,
    column_mapping: dict[str, str] | None,
    tag: str,
    ingest_id: int,
    start_index: int,
    expected_extra_columns: list[str],
) -> None:
    """
    Validate _process_cell_info_obs schema partitioning across mapping variants.

    :param small_csr_matrix: Small deterministic CSR matrix fixture
    :param column_mapping: Optional obs mapping dictionary
    :param tag: Tag value supplied to the serializer
    :param ingest_id: Ingest identifier used for the slice
    :param start_index: Starting index for the processed slice
    :param expected_extra_columns: Columns expected when mapping is provided
    """
    n_obs, _ = small_csr_matrix.shape
    obs = pd.DataFrame(
        data={
            "sample": ["A", "A", "B", "B"],
            "quality": [0.1, 0.9, 0.8, 0.2],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    adata = create_ingest_files.AnnData(X=small_csr_matrix, obs=obs)
    end_index = start_index + n_obs - 1

    df_schema, df_meta = create_ingest_files._process_cell_info_obs(
        adata=adata,
        tag=tag,
        ingest_id=ingest_id,
        start_index=start_index,
        end_index=end_index,
        column_mapping=column_mapping,
    )

    assert ingest_constants.OBS_NEXUS_ID in df_schema.columns
    assert ingest_constants.OBS_INGEST_ID in df_schema.columns
    assert ingest_constants.OBS_TAG in df_schema.columns
    assert ingest_constants.OBS_TOTAL_MRNA_UMIS in df_schema.columns
    assert df_schema[ingest_constants.OBS_TOTAL_MRNA_UMIS].tolist() == _expected_total_mrna_umis(adata.X)

    if expected_extra_columns:
        combined = pd.concat([df_schema, df_meta], axis=1)
        for column in expected_extra_columns:
            assert column in combined.columns


@pytest.mark.parametrize(
    "column_mapping, tag, ingest_id, start_index, expected_extra_columns",
    [
        pytest.param(None, "t1", 99, 1, [], id="without-mapping"),
        pytest.param(
            {"index": ingest_constants.VAR_FEATURE_INFO_ORIGINAL_ID, "gene": "gene_name"},
            "t2",
            77,
            10,
            [ingest_constants.VAR_FEATURE_INFO_ORIGINAL_ID, "gene_name"],
            id="with-mapping",
        ),
    ],
)
def test_process_feature_info_var_mapping_variants(
    column_mapping: dict[str, str] | None,
    tag: str,
    ingest_id: int,
    start_index: int,
    expected_extra_columns: list[str],
) -> None:
    """
    Validate _process_feature_info_var schema partitioning across mapping variants.

    :param column_mapping: Optional var mapping dictionary
    :param tag: Tag value supplied to the serializer
    :param ingest_id: Ingest identifier used for the slice
    :param start_index: Starting index for the processed slice
    :param expected_extra_columns: Columns expected when mapping is provided
    """
    df_var = pd.DataFrame(
        data={
            "gene": ["g1", "g2", "g3"],
            "chrom": ["1", "X", "MT"],
        }
    )
    end_index = start_index + len(df_var.index) - 1

    df_schema, df_meta = create_ingest_files._process_feature_info_var(
        df=df_var.copy(),
        tag=tag,
        ingest_id=ingest_id,
        start_index=start_index,
        end_index=end_index,
        column_mapping=column_mapping,
    )

    assert ingest_constants.VAR_NEXUS_ID in df_schema.columns
    assert ingest_constants.VAR_INGEST_ID in df_schema.columns
    assert ingest_constants.VAR_TAG in df_schema.columns

    if expected_extra_columns:
        combined = pd.concat([df_schema, df_meta], axis=1)
        for column in expected_extra_columns:
            assert column in combined.columns
