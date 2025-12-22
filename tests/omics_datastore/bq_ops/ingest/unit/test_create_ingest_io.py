import json
import pathlib

import anndata
import pandas as pd
import pytest
from fastavro import reader as avro_reader

from cellarium.nexus.omics_datastore.bq_ops import constants as ingest_constants
from cellarium.nexus.omics_datastore.bq_ops import exceptions as ingest_exceptions
from cellarium.nexus.omics_datastore.bq_ops.ingest import create_ingest_files


def _write_updated_h5ad_with_ids(src_path: pathlib.Path, dst_path: pathlib.Path) -> None:
    """
    Write a copy of the given .h5ad adding required id columns to obs/var.
    """
    ad = anndata.read_h5ad(filename=src_path)
    n_obs, n_vars = ad.shape
    ad.obs = ad.obs.copy()
    ad.var = ad.var.copy()
    ad.obs[ingest_constants.OBS_NEXUS_ID] = [f"c{i+1}" for i in range(n_obs)]
    ad.var[ingest_constants.VAR_NEXUS_ID] = [f"g{j+1}" for j in range(n_vars)]
    ad.write_h5ad(filename=dst_path)


def test_optimized_raw_matrix_read_coo_happy_path(small_anndata, tmp_path: pathlib.Path) -> None:
    _ad, src_path = small_anndata
    updated_path = tmp_path / "with_ids.h5ad"
    _write_updated_h5ad_with_ids(src_path=src_path, dst_path=updated_path)

    n_obs = anndata.read_h5ad(filename=updated_path).n_obs
    row_ids, col_ids, data = create_ingest_files.optimized_raw_matrix_read_coo(
        input_file_path=updated_path,
        row_offset=0,
        end=n_obs,
    )

    assert len(row_ids) == len(col_ids) == len(data)
    assert all(int(x) > 0 for x in data)
    assert all(str(r).startswith("c") for r in row_ids)
    assert all(str(c).startswith("g") for c in col_ids)


def test_optimized_raw_matrix_read_coo_missing_obs_id_raises(small_anndata, tmp_path: pathlib.Path) -> None:
    _ad, src_path = small_anndata
    dst = tmp_path / "missing_obs_id.h5ad"
    ad = anndata.read_h5ad(filename=src_path)
    n_obs, n_vars = ad.shape
    ad.obs = ad.obs.copy()
    ad.var = ad.var.copy()
    ad.var[ingest_constants.VAR_NEXUS_ID] = [f"g{j+1}" for j in range(n_vars)]
    ad.write_h5ad(filename=dst)

    with pytest.raises(KeyError):
        _ = create_ingest_files.optimized_raw_matrix_read_coo(
            input_file_path=dst,
            row_offset=0,
            end=n_obs,
        )


def test_optimized_raw_matrix_read_coo_missing_var_id_raises(small_anndata, tmp_path: pathlib.Path) -> None:
    _ad, src_path = small_anndata
    dst = tmp_path / "missing_var_id.h5ad"
    ad = anndata.read_h5ad(filename=src_path)
    n_obs, n_vars = ad.shape
    ad.obs = ad.obs.copy()
    ad.var = ad.var.copy()
    ad.obs[ingest_constants.OBS_NEXUS_ID] = [f"c{i+1}" for i in range(n_obs)]
    ad.write_h5ad(filename=dst)

    with pytest.raises(KeyError):
        _ = create_ingest_files.optimized_raw_matrix_read_coo(
            input_file_path=dst,
            row_offset=0,
            end=n_obs,
        )


def test_write_pandas_to_avro_info_file_round_trip(tmp_path: pathlib.Path, small_csr_matrix) -> None:
    n_obs, _ = small_csr_matrix.shape
    obs = pd.DataFrame(
        data={
            "sample": ["A", "A", "B", "B"],
            "quality": [0.1, 0.9, 0.8, 0.2],
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    ad = create_ingest_files.AnnData(X=small_csr_matrix, obs=obs)
    df_schema, df_meta = create_ingest_files._process_cell_info_obs(
        adata=ad,
        tag="tagA",
        ingest_id=1,
        start_index=1,
        end_index=n_obs,
        column_mapping={"index": ingest_constants.OBS_CELL_INFO_ORIGINAL_ID},
    )

    out_path = tmp_path / "cell-info.avro"
    create_ingest_files.write_pandas_to_avro_info_file(
        df=df_schema,
        df_metadata_extra=df_meta,
        pydantic_model_class=create_ingest_files.CellInfoBQAvroSchema,
        output_file_name=str(out_path),
    )

    with open(out_path, "rb") as f:
        records = list(avro_reader(f))

    assert len(records) == n_obs
    for rec in records:
        assert "metadata_extra" in rec
        meta = json.loads(rec["metadata_extra"]) if isinstance(rec["metadata_extra"], str) else rec["metadata_extra"]
        assert isinstance(meta, dict)


def test_dump_ingest_info_applies_size_limit_and_key_filtering(tmp_path: pathlib.Path, small_anndata) -> None:
    """
    Validate dump_ingest_info applies size-limit replacement for large values and supports key filtering.
    """
    ad, _ = small_anndata

    # Case 1: apply small size limit (so large_vec is replaced); keep all keys
    create_ingest_files.dump_ingest_info(
        ingest_id=42,
        adata_uns=ad.uns,
        output_dir=tmp_path,
        metadata_limit_size=64,  # intentionally tiny to trigger replacement
        uns_keys_to_keep=None,
    )

    ingest_path = tmp_path / ingest_constants.INGEST_INGEST_FILE_NAME
    with open(ingest_path, "rb") as f:
        records = list(avro_reader(f))
    assert len(records) == 1
    meta = (
        json.loads(records[0]["metadata_extra"])
        if isinstance(records[0]["metadata_extra"], str)
        else records[0]["metadata_extra"]
    )
    assert "large_vec" in meta
    assert isinstance(meta["large_vec"], str)
    assert "Value was removed due to large size" in meta["large_vec"]

    # Case 2: filter keys to only 'title'
    create_ingest_files.dump_ingest_info(
        ingest_id=42,
        adata_uns=ad.uns,
        output_dir=tmp_path,
        metadata_limit_size=10_000_000,
        uns_keys_to_keep=["title"],
    )
    with open(ingest_path, "rb") as f:
        records2 = list(avro_reader(f))
    # File is opened in append mode; read the last record to validate the most recent write
    last = records2[-1]
    meta2 = json.loads(last["metadata_extra"]) if isinstance(last["metadata_extra"], str) else last["metadata_extra"]
    assert set(meta2.keys()) == {"title"}


def test_dump_core_matrix_batch_writes_expected_csv(tmp_path: pathlib.Path, small_anndata) -> None:
    """
    Validate that dump_core_matrix_batch writes expected COO triples to CSV for a small slice.
    """
    _ad, src_path = small_anndata
    updated_path = tmp_path / "with_ids_for_batch.h5ad"
    _write_updated_h5ad_with_ids(src_path=src_path, dst_path=updated_path)

    ad2 = anndata.read_h5ad(filename=updated_path)
    n_obs, _ = ad2.shape
    out_csv = tmp_path / "batch.csv"

    create_ingest_files.dump_core_matrix_batch(
        input_file_path=updated_path,
        row_offset=0,
        end=n_obs,
        batch_output_file_name=str(out_csv),
        batch_num=0,
    )

    # Build expected set of lines from the file content itself (ensuring non-zero entries only)
    coord = ad2.X[:].tocoo()
    row_ids = ad2.obs[ingest_constants.OBS_NEXUS_ID].iloc[coord.row + 0].tolist()
    col_ids = ad2.var[ingest_constants.VAR_NEXUS_ID].iloc[coord.col].tolist()
    data_vals = coord.data.tolist()
    expected_lines = [f"{r},{c},{int(v)}\n" for r, c, v in zip(row_ids, col_ids, data_vals) if int(v) != 0]

    with open(out_csv, "r") as fh:
        written_lines = fh.readlines()

    assert sorted(written_lines) == sorted(expected_lines)


def test_dump_core_matrix_in_parallel_creates_batches_and_propagates_errors(
    tmp_path: pathlib.Path, small_anndata, patched_process_pool_executor, monkeypatch: pytest.MonkeyPatch
) -> None:
    """
    Validate batch file creation with patched executor and error propagation to DataIngestError.
    """
    _ad, src_path = small_anndata
    updated_path = tmp_path / "with_ids_parallel.h5ad"
    _write_updated_h5ad_with_ids(src_path=src_path, dst_path=updated_path)

    n_obs = anndata.read_h5ad(filename=updated_path).n_obs

    # Success path: should create one batch file for small n_obs
    create_ingest_files.dump_core_matrix_in_parallel(
        input_file_path=updated_path,
        total_cells=n_obs,
        output_dir=tmp_path,
    )
    created = list(tmp_path.glob(ingest_constants.INGEST_RAW_COUNTS_FILE_PATTERN))
    assert len(created) >= 1

    # Error propagation: make batch function raise and expect DataIngestError
    def _raise(*args, **kwargs):  # noqa: ANN001, D401
        raise RuntimeError("boom")

    monkeypatch.setattr(create_ingest_files, "dump_core_matrix_batch", _raise)
    with pytest.raises(ingest_exceptions.DataIngestError):
        create_ingest_files.dump_core_matrix_in_parallel(
            input_file_path=updated_path,
            total_cells=n_obs,
            output_dir=tmp_path,
        )
