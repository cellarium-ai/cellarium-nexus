import json
import pathlib

import anndata
import pytest
from fastavro import reader as avro_reader

from cellarium.nexus.omics_datastore.bq_ops import constants as ingest_constants
from cellarium.nexus.omics_datastore.bq_ops import exceptions as ingest_exceptions
from cellarium.nexus.omics_datastore.bq_ops.ingest import create_ingest_files


def test_create_ingest_files_with_mapping_and_filtering(
    tmp_path: pathlib.Path,
    small_anndata: tuple[anndata.AnnData, pathlib.Path],
    patched_process_pool_executor: None,  # ensures thread-based parallelism for deterministic tests
) -> None:
    """
    Exercise end-to-end creation of ingest artifacts with column mapping and uns key filtering.
    """
    ad, ad_path = small_anndata
    n_obs, n_vars = ad.shape

    column_mapping = {
        "obs_mapping": {
            "index": ingest_constants.OBS_CELL_INFO_ORIGINAL_ID,
            "sample": "sample_name",
        },
        # Map required FeatureInfo string fields as well
        # index -> ensemble_id, gene -> symbol, chrom -> reference
        "var_mapping": {
            "index": "ensemble_id",
            "gene": "symbol",
            "chrom": "reference",
        },
    }

    result = create_ingest_files.create_ingest_files(
        adata_file_path=ad_path,
        tag="e2e_tag",
        cell_info_start_index=1,
        cell_info_end_index=n_obs,
        feature_info_start_index=1,
        feature_info_end_index=n_vars,
        ingest_id=101,
        output_dir=tmp_path,
        column_mapping=column_mapping,
        uns_keys_to_keep=["title"],
    )

    # Validate result structure
    assert result["output_dir"] == str(tmp_path)
    assert isinstance(result["adata_uns_clean"], dict)

    # Validate files exist
    ingest_info_path = tmp_path / ingest_constants.INGEST_INGEST_FILE_NAME
    cell_info_path = tmp_path / ingest_constants.INGEST_CELL_INFO_FILE_NAME
    feature_info_path = tmp_path / ingest_constants.INGEST_FEATURE_INFO_FILE_NAME
    raw_counts_files = list(tmp_path.glob(ingest_constants.INGEST_RAW_COUNTS_FILE_PATTERN))

    assert ingest_info_path.exists()
    assert cell_info_path.exists()
    assert feature_info_path.exists()
    assert len(raw_counts_files) >= 1

    # Validate ingest-info content (check last record for filtering)
    with open(ingest_info_path, "rb") as f:
        records = list(avro_reader(f))
    assert len(records) >= 1
    last = records[-1]
    assert last["id"] == 101
    meta = json.loads(last["metadata_extra"]) if isinstance(last["metadata_extra"], str) else last["metadata_extra"]
    assert set(meta.keys()) == {"title"}

    # Validate basic structure of cell-info and feature-info
    with open(cell_info_path, "rb") as f:
        cell_records = list(avro_reader(f))
    with open(feature_info_path, "rb") as f:
        feature_records = list(avro_reader(f))

    assert len(cell_records) == n_obs
    assert len(feature_records) == n_vars

    # Spot-check tag and ids present
    assert all("id" in r and "ingest_id" in r and "tag" in r for r in cell_records)
    assert all("id" in r and "ingest_id" in r and "tag" in r for r in feature_records)

    # Spot-check one raw-counts CSV format
    sample_csv = raw_counts_files[0]
    with open(sample_csv, "r") as fh:
        line = fh.readline().strip()
    # If there are no non-zero entries in first file, skip check
    if line:
        parts = line.split(",")
        assert len(parts) == 3
        assert parts[2].isdigit()


def test_create_ingest_files_tag_empty_string(
    tmp_path: pathlib.Path,
    small_anndata: tuple[anndata.AnnData, pathlib.Path],
    patched_process_pool_executor: None,
) -> None:
    """
    Exercise end-to-end creation with tag=None to ensure Avro tag fields can be null.
    """
    ad, ad_path = small_anndata
    n_obs, n_vars = ad.shape

    # Provide mappings to satisfy validators
    column_mapping = {
        "obs_mapping": {
            "index": ingest_constants.OBS_CELL_INFO_ORIGINAL_ID,
        },
        "var_mapping": {
            "index": "ensemble_id",
            "gene": "symbol",
            "chrom": "reference",
        },
    }

    result = create_ingest_files.create_ingest_files(
        adata_file_path=ad_path,
        tag="",
        cell_info_start_index=1,
        cell_info_end_index=n_obs,
        feature_info_start_index=1,
        feature_info_end_index=n_vars,
        ingest_id=202,
        output_dir=tmp_path,
        column_mapping=column_mapping,
    )

    ingest_info_path = tmp_path / ingest_constants.INGEST_INGEST_FILE_NAME
    cell_info_path = tmp_path / ingest_constants.INGEST_CELL_INFO_FILE_NAME
    feature_info_path = tmp_path / ingest_constants.INGEST_FEATURE_INFO_FILE_NAME
    raw_counts_files = list(tmp_path.glob(ingest_constants.INGEST_RAW_COUNTS_FILE_PATTERN))

    assert ingest_info_path.exists()
    assert cell_info_path.exists()
    assert feature_info_path.exists()
    assert len(raw_counts_files) >= 1

    with open(cell_info_path, "rb") as f:
        cell_records = list(avro_reader(f))
    with open(feature_info_path, "rb") as f:
        feature_records = list(avro_reader(f))

    assert len(cell_records) == n_obs
    assert len(feature_records) == n_vars

    # Tag stored as empty string
    assert all((r.get("tag") == "" or r.get("tag") is None) for r in cell_records)
    assert all((r.get("tag") == "" or r.get("tag") is None) for r in feature_records)


def test_create_ingest_files_raises_when_required_string_columns_missing(
    tmp_path: pathlib.Path,
    small_anndata: tuple[anndata.AnnData, pathlib.Path],
    patched_process_pool_executor: None,
) -> None:
    """
    Expect DataIngestError when required string columns are missing (e.g., obs original_id, var ensemble_id/reference).
    """
    _ad, ad_path = small_anndata
    # Do not provide column_mapping so required fields remain missing/NA after schema ensure
    with pytest.raises(ingest_exceptions.DataValidationError) as exc:
        _ = create_ingest_files.create_ingest_files(
            adata_file_path=ad_path,
            tag="e2e",
            cell_info_start_index=1,
            cell_info_end_index=4,
            feature_info_start_index=1,
            feature_info_end_index=3,
            ingest_id=303,
            output_dir=tmp_path,
        )
    assert "Missing required" in str(exc.value) or "Invalid values for required" in str(exc.value)
