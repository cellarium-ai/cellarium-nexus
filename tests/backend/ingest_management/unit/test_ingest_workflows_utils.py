from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from django.conf import settings as dj_settings

from cellarium.nexus.backend.cell_management import models as cm_models
from cellarium.nexus.backend.ingest_management.utils import workflows_utils


@pytest.mark.usefixtures("vertex_ai_pipeline_stub")
@pytest.mark.django_db
def test_submit_ingest_pipeline_creates_configs_and_calls_submit(
    default_dataset: cm_models.OmicsDataset, monkeypatch: pytest.MonkeyPatch, vertex_ai_pipeline_stub
) -> None:
    """
    Ensure submit_ingest_pipeline builds one CreateIngestFilesConfig per input row,
    dumps both create+ingest configs, and calls submit_pipeline with their paths.
    """
    # Arrange input dataframe (two files, with tags)
    df = pd.DataFrame(
        [
            {"gcs_file_path": "gs://bucket/data/a1.h5ad", "tag": "t1", "ingest_id": 10, "ingest_file_id": 101},
            {"gcs_file_path": "gs://bucket/data/b2.h5ad", "tag": "t2", "ingest_id": 10, "ingest_file_id": 102},
        ]
    )

    # Capture configs passed to dump
    captured: dict[str, Any] = {"create": None, "ingest": None}

    def _dump_configs(configs, bucket_path, workers=8):
        # Determine type by first element's class name
        if not configs:
            return []
        model_name = configs[0].__class__.__name__
        if model_name.lower().startswith("createingestfilesconfig"):
            captured["create"] = configs
            return [f"gs://cfg/create_{i}.yaml" for i in range(len(configs))]
        else:
            captured["ingest"] = configs
            return ["gs://cfg/ingest.yaml"]

    monkeypatch.setattr(
        workflows_utils.utils.workflows_configs,
        "dump_configs_to_bucket",
        _dump_configs,
        raising=True,
    )

    # Act
    url = workflows_utils.submit_ingest_pipeline(
        df_ingest_file_info=df,
        omics_dataset=default_dataset,
        column_mapping={"donor_id": "donor_id"},
        validation_methods=["validate_shapes"],
    )

    # Assert pipeline URL from stub
    assert isinstance(url, str) and url == vertex_ai_pipeline_stub.default_url

    # Assert we created one CreateIngestFilesConfig per row
    assert captured["create"] is not None and len(captured["create"]) == 2
    for cfg, row in zip(captured["create"], df.to_dict(orient="records")):
        assert cfg.bigquery_dataset == default_dataset.name
        assert cfg.input_file_path == row["gcs_file_path"]
        assert cfg.tag == row["tag"]
        assert cfg.ingest_id == row["ingest_id"]
        assert cfg.ingest_file_id == row["ingest_file_id"]
        assert cfg.column_mapping == {"donor_id": "donor_id"}
        assert cfg.validation_methods == ["validate_shapes"]
        # Unique stage dirs, ending with filename stem prefix
        stem = row["gcs_file_path"].split("/")[-1].split(".")[0][:10]
        assert cfg.bucket_stage_dir.endswith(stem)
    # Ingest config should aggregate stage dirs from create configs
    assert captured["ingest"] is not None and len(captured["ingest"]) == 1
    ingest_cfg = captured["ingest"][0]
    assert len(ingest_cfg.bucket_stage_dirs) == 2
    assert all(sd in ingest_cfg.bucket_stage_dirs for sd in [c.bucket_stage_dir for c in captured["create"]])


@pytest.mark.usefixtures("vertex_ai_pipeline_stub")
def test_submit_validation_pipeline_batches_and_calls_submit(
    monkeypatch: pytest.MonkeyPatch, vertex_ai_pipeline_stub
) -> None:
    """
    Ensure submit_validation_pipeline batches adata paths by MAX_ADATA_FILES_PER_VALIDATION_BATCH,
    dumps all ValidationConfig batches, and calls submit_pipeline with the resulting paths.
    """
    # Force small batch size to exercise batching
    monkeypatch.setattr(dj_settings, "MAX_ADATA_FILES_PER_VALIDATION_BATCH", 3, raising=False)

    adata_paths = [f"gs://bucket/d{i}.h5ad" for i in range(7)]  # -> 3,3,1 batches
    validation_methods = ["validate_shapes", "validate_schema"]

    captured: dict[str, Any] = {"validation": None}

    def _dump_configs(configs, bucket_path, workers=8):
        captured["validation"] = configs
        # Return deterministic paths per config
        return [f"gs://cfg/validation_{i}.yaml" for i in range(len(configs))]

    monkeypatch.setattr(
        workflows_utils.utils.workflows_configs,
        "dump_configs_to_bucket",
        _dump_configs,
        raising=True,
    )

    url = workflows_utils.submit_validation_pipeline(
        adata_gcs_paths=adata_paths,
        validation_report_id=123,
        validation_methods=validation_methods,
    )

    assert isinstance(url, str) and url == vertex_ai_pipeline_stub.default_url

    # We expect 3 configs (3,3,1)
    assert captured["validation"] is not None
    configs = captured["validation"]
    assert len(configs) == 3
    # Check payloads in each config
    assert configs[0].adata_gcs_paths == adata_paths[:3]
    assert configs[1].adata_gcs_paths == adata_paths[3:6]
    assert configs[2].adata_gcs_paths == adata_paths[6:]
    for cfg in configs:
        assert cfg.validation_report_id == 123
        assert cfg.validation_methods == validation_methods
        assert cfg.max_bytes_valid_per_file == dj_settings.INGEST_INPUT_FILE_MAX_SIZE
