from __future__ import annotations

import sys
import types
from typing import Any

import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.utils import workflows_utils


def _make_feature_schema(*, symbols: list[tuple[str, str]] = None) -> models.FeatureSchema:
    """
    Create a FeatureSchema with associated FeatureInfo records for testing.

    :param symbols: List of (ensemble_id, symbol) pairs

    :return: Persisted FeatureSchema instance with features populated
    """
    fs = models.FeatureSchema.objects.create(name="schema-test")
    for ens, sym in symbols or [("ENSG000001", "TP53"), ("ENSG000002", "EGFR")]:
        models.FeatureInfo.objects.create(
            ensemble_id=ens,
            symbol=sym,
            feature_schema=fs,
        )
    return fs


def test_get_total_cell_in_bq_number_counts(
    monkeypatch: pytest.MonkeyPatch, default_dataset: models.BigQueryDataset
) -> None:
    """
    Assert BigQueryDataOperator.count_cells is invoked and its value returned.
    """

    class _Op:
        def __init__(self, *args, **kwargs) -> None:
            pass

        def count_cells(self, *args, **kwargs) -> int:
            return 1234

    monkeypatch.setattr(workflows_utils, "BigQueryDataOperator", _Op, raising=True)

    out = workflows_utils.get_total_cell_in_bq_number(bigquery_dataset=default_dataset, filters={"a": 1})
    assert out == 1234


def test_compose_extract_curriculum_configs_invalid_bin_size(default_dataset: models.BigQueryDataset) -> None:
    """
    Raise ValueError when extract_bin_size <= 0.
    """
    fs = _make_feature_schema()

    with pytest.raises(ValueError):
        workflows_utils.compose_extract_curriculum_configs(
            name="x",
            creator_id=1,
            feature_schema=fs,
            bigquery_dataset=default_dataset,
            extract_bin_size=0,
            categorical_column_count_limit=2000,
            obs_columns=["organism"],
        )


def test_compose_extract_curriculum_configs_zero_cells(
    monkeypatch: pytest.MonkeyPatch, default_dataset: models.BigQueryDataset
) -> None:
    """
    Raise ZeroCellsReturnedError when total cell count is zero.
    """
    fs = _make_feature_schema()
    monkeypatch.setattr(workflows_utils, "get_total_cell_in_bq_number", lambda **_: 0, raising=True)

    with pytest.raises(workflows_utils.exceptions.ZeroCellsReturnedError):
        workflows_utils.compose_extract_curriculum_configs(
            name="x",
            creator_id=1,
            feature_schema=fs,
            bigquery_dataset=default_dataset,
            extract_bin_size=100,
            categorical_column_count_limit=2000,
            obs_columns=["organism"],
        )


def test_compose_extract_curriculum_configs_happy_path(
    monkeypatch: pytest.MonkeyPatch, default_dataset: models.BigQueryDataset
) -> None:
    """
    Compose configs with correct number of bins, features, and merged obs columns.
    """
    fs = _make_feature_schema(symbols=[("E1", "S1"), ("E2", "S2"), ("E3", "S3")])
    # total_cells=5, bin_size=2 => 3 bins
    monkeypatch.setattr(workflows_utils, "get_total_cell_in_bq_number", lambda **_: 5, raising=True)

    prepare_cfg, extract_cfgs = workflows_utils.compose_extract_curriculum_configs(
        name="extract-1",
        creator_id=7,
        feature_schema=fs,
        bigquery_dataset=default_dataset,
        extract_bin_size=2,
        categorical_column_count_limit=2000,
        obs_columns=["organism"],
        metadata_extra_columns=["extra1"],
    )

    # Features mapped 1:1 from schema
    assert len(prepare_cfg.features) == fs.features.count()
    # With default BINS_PER_WORKER (32), all bins in one worker config
    assert len(extract_cfgs) == 1
    cfg = extract_cfgs[0]
    assert cfg.bins == [0, 1, 2]
    # metadata_extra merged into obs columns for extract
    assert "extra1" in cfg.obs_columns


def test_compose_and_dump_configs_delegates_and_returns_paths(
    monkeypatch: pytest.MonkeyPatch, default_dataset: models.BigQueryDataset
) -> None:
    """
    Ensure compose_and_dump_configs dumps both prepare and extract configs and returns their paths.
    """
    fs = _make_feature_schema()
    monkeypatch.setattr(workflows_utils, "get_total_cell_in_bq_number", lambda **_: 10, raising=True)

    captured: dict[str, Any] = {"calls": []}

    def _dump(configs, bucket_path, workers=8):
        # record number of configs by type
        captured["calls"].append((type(configs[0]).__name__, len(configs)))
        # return predictable paths
        if type(configs[0]).__name__.lower().startswith("bqopsprepareextract"):
            return ["gs://bucket/prepare.yaml"]
        else:
            return [f"gs://bucket/extract_{i}.yaml" for i in range(len(configs))]

    monkeypatch.setattr(
        workflows_utils.workflows_configs,
        "dump_configs_to_bucket",
        _dump,
        raising=True,
    )

    prepare_path, extract_paths = workflows_utils.compose_and_dump_configs(
        feature_schema=fs,
        creator_id=42,
        bigquery_dataset=default_dataset,
        name="ex",
        extract_bin_size=5,
        categorical_column_count_limit=2000,
        obs_columns=["organism"],
    )

    assert prepare_path == "gs://bucket/prepare.yaml"
    assert all(p.startswith("gs://bucket/extract_") for p in extract_paths)
    # Ensure we dumped 1 prepare and N extract configs
    names = [n for n, _ in captured["calls"]]
    assert "BQOpsPrepareExtract" in names
    assert any(n == "BQOpsExtract" for n in names)
    assert any(c > 0 for n, c in captured["calls"] if n == "BQOpsExtract")


@pytest.mark.usefixtures("vertex_ai_pipeline_stub")
def test_submit_extract_pipeline_calls_submit(
    monkeypatch: pytest.MonkeyPatch,
    default_dataset: models.BigQueryDataset,
    vertex_ai_pipeline_stub,
) -> None:
    """
    Ensure submit_extract_pipeline calls submit_pipeline with expected kwargs.
    """
    fs = _make_feature_schema()

    # Short-circuit heavy work; assume compose returns paths
    monkeypatch.setattr(
        workflows_utils,
        "compose_and_dump_configs",
        lambda **kwargs: ("gs://bucket/prepare.yaml", ["gs://bucket/extract_0.yaml", "gs://bucket/extract_1.yaml"]),
        raising=True,
    )

    # Provide a dummy kubeflow pipelines module to satisfy import inside submit_extract_pipeline
    dummy_pipelines_mod = types.SimpleNamespace(extract_data_pipeline=object())
    monkeypatch.setitem(
        sys.modules,
        "cellarium.nexus.workflows.kubeflow.pipelines",
        dummy_pipelines_mod,
    )

    url = workflows_utils.submit_extract_pipeline(
        feature_schema=fs,
        creator_id=1,
        bigquery_dataset=default_dataset,
        name="e",
        extract_bin_size=10,
        categorical_column_count_limit=2000,
        obs_columns=["organism"],
    )

    # VertexAIPipelineStub returns a default URL and records the last call
    assert url == vertex_ai_pipeline_stub.default_url

    # Verify submit_pipeline call payload
    call = vertex_ai_pipeline_stub.calls[-1]
    kwargs = call["kwargs"]
    assert kwargs["pipeline_component"] == dummy_pipelines_mod.extract_data_pipeline
    assert kwargs["pipeline_kwargs"] == {
        "prepare_extract_config": "gs://bucket/prepare.yaml",
        "extract_configs": [
            "gs://bucket/extract_0.yaml",
            "gs://bucket/extract_1.yaml",
        ],
    }
