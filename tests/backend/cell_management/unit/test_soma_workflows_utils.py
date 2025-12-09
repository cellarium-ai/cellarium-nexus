"""
Unit tests for SOMA workflows utilities (function-based).
"""

from unittest import mock

import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.utils import exceptions, soma_workflows_utils
from cellarium.nexus.shared.schemas import component_configs


@pytest.fixture
def soma_dataset() -> models.OmicsDataset:
    dataset = mock.Mock(spec=models.OmicsDataset)
    dataset.name = "test_soma_dataset"
    dataset.uri = "gs://bucket/soma_experiment"
    dataset.backend = models.OmicsDatasetBackend.TILEDB_SOMA
    return dataset


def _make_settings(**overrides) -> mock.Mock:
    settings = mock.Mock()
    settings.BACKEND_PIPELINE_DIR = "backend"
    settings.PIPELINE_DATA_EXTRACTS_DIR = "extracts"
    settings.SITE_URL = "https://nexus.example.com"
    settings.BUCKET_NAME_PRIVATE = "private-bucket"
    settings.PIPELINE_CONFIGS_DIR = "configs"
    settings.GCP_PROJECT_ID = "test-project"
    settings.PIPELINE_SERVICE_ACCOUNT = "sa@test.iam.gserviceaccount.com"
    settings.PIPELINE_ROOT_PATH = "gs://bucket/pipeline-root"
    settings.GCP_APPLICATION_BILLING_LABEL = "nexus"
    settings.TILEDB_SOMA_RANGES_PER_WORKER = 32
    for key, value in overrides.items():
        setattr(settings, key, value)
    return settings


@pytest.fixture(autouse=True)
def _stub_tiledb_operator(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Prevent expensive TileDB interactions during unit tests.

    Provide a lightweight operator with a no-op count_cells to keep tests fast.
    """

    class _Op:
        def count_cells(self, *, filter_statements: dict | None = None) -> int:
            return 0

    monkeypatch.setattr(soma_workflows_utils, "TileDBSOMADataOperator", lambda experiment_uri=None: _Op(), raising=True)


@pytest.fixture(autouse=True)
def _stub_coordinator(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Stub SomaDataOpsCoordinator to avoid real I/O or TileDB calls during unit tests.
    """

    class _Coordinator:
        def __init__(self, *args, **kwargs) -> None:
            self.prepare_soma_extract = mock.Mock()
            self.run_soma_extract = mock.Mock()
            self.mark_soma_curriculum_as_finished = mock.Mock()

    monkeypatch.setattr(soma_workflows_utils, "SomaDataOpsCoordinator", lambda **kwargs: _Coordinator(), raising=True)


def test_get_total_cells_in_soma_counts(monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset) -> None:
    filters = {"cell_type__in": ["T cell", "B cell"]}

    class _Op:
        def count_cells(self, *, filter_statements: dict | None = None) -> int:
            return 1000

    monkeypatch.setattr(soma_workflows_utils, "TileDBSOMADataOperator", lambda experiment_uri: _Op(), raising=True)

    result = soma_workflows_utils.get_total_cells_in_soma(omics_dataset=soma_dataset, filters=filters)

    assert result == 1000


def test_get_total_cells_in_soma_raises_without_uri(soma_dataset: models.OmicsDataset) -> None:
    soma_dataset.uri = None

    with pytest.raises(ValueError, match="has no URI configured"):
        soma_workflows_utils.get_total_cells_in_soma(omics_dataset=soma_dataset)


def test_compose_soma_extract_configs_invalid_range_size(soma_dataset: models.OmicsDataset) -> None:
    with pytest.raises(ValueError, match="Range size must be greater than 0"):
        soma_workflows_utils.compose_soma_extract_configs(
            name="test_extract",
            creator_id=1,
            omics_dataset=soma_dataset,
            range_size=0,
            output_chunk_size=1000,
        )


def test_compose_soma_extract_configs_no_uri(soma_dataset: models.OmicsDataset) -> None:
    soma_dataset.uri = None

    with pytest.raises(ValueError, match="has no URI configured"):
        soma_workflows_utils.compose_soma_extract_configs(
            name="test_extract",
            creator_id=1,
            omics_dataset=soma_dataset,
            range_size=1000,
            output_chunk_size=1000,
        )


def test_compose_soma_extract_configs_zero_cells(
    monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset
) -> None:
    monkeypatch.setattr(soma_workflows_utils, "get_total_cells_in_soma", lambda **kwargs: 0, raising=True)

    with pytest.raises(exceptions.ZeroCellsReturnedError):
        soma_workflows_utils.compose_soma_extract_configs(
            name="test_extract",
            creator_id=1,
            omics_dataset=soma_dataset,
            range_size=1000,
            output_chunk_size=1000,
        )


def test_compose_soma_extract_configs_happy_path(
    monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset
) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)
    monkeypatch.setattr(soma_workflows_utils, "get_total_cells_in_soma", lambda **kwargs: 5000, raising=True)

    extract_configs = soma_workflows_utils.compose_soma_extract_configs(
        name="test_extract",
        creator_id=42,
        omics_dataset=soma_dataset,
        range_size=1000,
        output_chunk_size=1000,
        filters={"cell_type__eq": "T cell"},
    )

    assert len(extract_configs) == 1
    cfg = extract_configs[0]
    assert isinstance(cfg, component_configs.SomaOpsExtract)
    assert cfg.extract_name == "test_extract"
    assert cfg.plan_path == "backend/extracts/test_extract/soma_extract_plan.json"
    assert cfg.range_indices_slice == [0, 1, 2, 3, 4]


def test_compose_soma_extract_configs_multiple_batches(
    monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset
) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)
    monkeypatch.setattr(soma_workflows_utils, "get_total_cells_in_soma", lambda **kwargs: 100_000, raising=True)

    extract_configs = soma_workflows_utils.compose_soma_extract_configs(
        name="test_extract",
        creator_id=1,
        omics_dataset=soma_dataset,
        range_size=1000,  # 100 ranges
        output_chunk_size=1000,
    )

    assert len(extract_configs) == 4
    assert extract_configs[0].range_indices_slice == list(range(0, 32))
    assert extract_configs[1].range_indices_slice == list(range(32, 64))
    assert extract_configs[2].range_indices_slice == list(range(64, 96))
    assert extract_configs[3].range_indices_slice == list(range(96, 100))


def test_compose_and_dump_soma_configs(monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)

    mock_extract_configs = [
        mock.Mock(
            experiment_uri="gs://exp",
            nexus_backend_api_url="https://site",
            bucket_name="private-bucket",
            plan_path="plan1",
            extract_name="test_extract",
        ),
        mock.Mock(
            experiment_uri="gs://exp",
            nexus_backend_api_url="https://site",
            bucket_name="private-bucket",
            plan_path="plan2",
            extract_name="test_extract",
        ),
    ]
    monkeypatch.setattr(
        soma_workflows_utils,
        "compose_soma_extract_configs",
        lambda **kwargs: mock_extract_configs,
        raising=True,
    )

    coordinator = mock.Mock()
    monkeypatch.setattr(soma_workflows_utils, "SomaDataOpsCoordinator", lambda **kwargs: coordinator, raising=True)

    monkeypatch.setattr(
        soma_workflows_utils.workflows_configs,
        "dump_configs_to_bucket",
        lambda configs, bucket_path: [f"gs://bucket/extract_{i}.yaml" for i, _ in enumerate(configs)],
        raising=True,
    )

    paths = soma_workflows_utils.compose_and_dump_soma_configs(
        name="test_extract",
        creator_id=1,
        omics_dataset=soma_dataset,
        range_size=1000,
        output_chunk_size=1000,
        prepare=False,
    )

    assert paths == ["gs://bucket/extract_0.yaml", "gs://bucket/extract_1.yaml"]
    coordinator.prepare_soma_extract.assert_not_called()


def test_submit_soma_extract_pipeline(monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)

    monkeypatch.setattr(
        soma_workflows_utils,
        "compose_and_dump_soma_configs",
        lambda **kwargs: ["gs://bucket/extract_1.yaml"],
        raising=True,
    )

    submit_stub = mock.Mock(return_value="https://console.cloud.google.com/vertex-ai/pipelines/123")
    monkeypatch.setattr(
        "cellarium.nexus.workflows.kubeflow.utils.job.submit_pipeline",
        submit_stub,
        raising=True,
    )

    url = soma_workflows_utils.submit_soma_extract_pipeline(
        name="test_extract",
        creator_id=1,
        omics_dataset=soma_dataset,
        range_size=1000,
        output_chunk_size=1000,
    )

    submit_stub.assert_called_once()
    kwargs = submit_stub.call_args.kwargs
    assert kwargs["display_name"] == "Nexus SOMA Extract - test_extract"
    assert kwargs["gcp_project"] == "test-project"
    assert kwargs["pipeline_kwargs"]["extract_configs"] == ["gs://bucket/extract_1.yaml"]
    assert kwargs["pipeline_kwargs"]["mark_finished_config"] == "gs://bucket/extract_1.yaml"
    assert url == "https://console.cloud.google.com/vertex-ai/pipelines/123"
