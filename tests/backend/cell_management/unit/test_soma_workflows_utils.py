"""
Unit tests for SOMA workflows utilities (function-based).
"""

from unittest import mock

import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.utils import soma_workflows_utils
from cellarium.nexus.shared.schemas import component_configs
from cellarium.nexus.shared.schemas.omics_datastore import IdContiguousRange, RandomizedCurriculumMetadata


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


def test_compose_soma_randomized_extract_configs_no_uri(soma_dataset: models.OmicsDataset) -> None:
    soma_dataset.uri = None
    curriculum_metadata = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma_experiment",
        value_filter="",
        total_cells=5000,
        id_ranges=[IdContiguousRange(start=0, end=1000)],
        range_size=1000,
        num_ranges=1,
        num_bins=5,
        extract_bin_size=1000,
        last_bin_size=1000,
        extract_bin_indexes=[0, 1, 2, 3, 4],
    )

    with pytest.raises(ValueError, match="has no URI configured"):
        soma_workflows_utils.compose_soma_randomized_extract_configs(
            name="test_extract",
            omics_dataset=soma_dataset,
            curriculum_metadata=curriculum_metadata,
            extract_metadata_path="path/to/metadata.json",
            extract_bucket_path="path/to/extract",
        )


def test_compose_soma_randomized_extract_configs_happy_path(
    monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset
) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)

    # 5 ranges, 32 ranges per worker -> 1 worker
    curriculum_metadata = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma_experiment",
        value_filter="cell_type == 'T cell'",
        total_cells=5000,
        id_ranges=[IdContiguousRange(start=i * 1000, end=(i + 1) * 1000) for i in range(5)],
        range_size=1000,
        num_ranges=5,
        num_bins=5,
        extract_bin_size=1000,
        last_bin_size=1000,
        extract_bin_indexes=[0, 1, 2, 3, 4],
    )

    extract_configs = soma_workflows_utils.compose_soma_randomized_extract_configs(
        name="test_extract",
        omics_dataset=soma_dataset,
        curriculum_metadata=curriculum_metadata,
        extract_metadata_path="backend/extracts/test_extract/extract_metadata.json",
        extract_bucket_path="backend/extracts/test_extract",
    )

    assert len(extract_configs) == 1
    cfg = extract_configs[0]
    assert isinstance(cfg, component_configs.SomaOpsExtract)
    assert cfg.extract_name == "test_extract"
    assert cfg.extract_metadata_path == "backend/extracts/test_extract/extract_metadata.json"
    assert cfg.partition_index == 0
    assert cfg.max_ranges_per_partition == 32  # From settings.TILEDB_SOMA_RANGES_PER_WORKER


def test_compose_soma_randomized_extract_configs_multiple_workers(
    monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset
) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)

    # 100 ranges, 32 ranges per worker -> 4 workers
    curriculum_metadata = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma_experiment",
        value_filter="",
        total_cells=100_000,
        id_ranges=[IdContiguousRange(start=i * 1000, end=(i + 1) * 1000) for i in range(100)],
        range_size=1000,
        num_ranges=100,
        num_bins=100,
        extract_bin_size=1000,
        last_bin_size=1000,
        extract_bin_indexes=list(range(100)),
    )

    extract_configs = soma_workflows_utils.compose_soma_randomized_extract_configs(
        name="test_extract",
        omics_dataset=soma_dataset,
        curriculum_metadata=curriculum_metadata,
        extract_metadata_path="path/to/metadata.json",
        extract_bucket_path="path/to/extract",
    )

    assert len(extract_configs) == 4
    for i, cfg in enumerate(extract_configs):
        assert cfg.partition_index == i
        assert cfg.max_ranges_per_partition == 32  # From settings.TILEDB_SOMA_RANGES_PER_WORKER


def test_compose_and_dump_soma_configs(monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset) -> None:
    settings = _make_settings()
    monkeypatch.setattr(soma_workflows_utils, "settings", settings, raising=False)

    mock_extract_configs = [
        mock.Mock(
            experiment_uri="gs://exp",
            nexus_backend_api_url="https://site",
            bucket_name="private-bucket",
            extract_metadata_path="metadata1",
            extract_name="test_extract",
        ),
        mock.Mock(
            experiment_uri="gs://exp",
            nexus_backend_api_url="https://site",
            bucket_name="private-bucket",
            extract_metadata_path="metadata2",
            extract_name="test_extract",
        ),
    ]
    monkeypatch.setattr(
        soma_workflows_utils,
        "compose_soma_randomized_extract_configs",
        lambda **kwargs: mock_extract_configs,
        raising=True,
    )

    curriculum_metadata = RandomizedCurriculumMetadata(
        experiment_uri="gs://bucket/soma_experiment",
        value_filter="",
        total_cells=5000,
        id_ranges=[IdContiguousRange(start=0, end=1000)],
        range_size=1000,
        num_ranges=1,
        num_bins=5,
        extract_bin_size=1000,
        last_bin_size=1000,
        extract_bin_indexes=[0, 1, 2, 3, 4],
    )
    coordinator = mock.Mock()
    coordinator.prepare_soma_extract.return_value = curriculum_metadata
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
        extract_bin_size=1000,
    )

    assert paths == ["gs://bucket/extract_0.yaml", "gs://bucket/extract_1.yaml"]
    coordinator.prepare_soma_extract.assert_called_once()


def test_submit_soma_randomized_extract_pipeline(
    monkeypatch: pytest.MonkeyPatch, soma_dataset: models.OmicsDataset
) -> None:
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

    url = soma_workflows_utils.submit_soma_randomized_extract_pipeline(
        name="test_extract",
        creator_id=1,
        omics_dataset=soma_dataset,
        range_size=1000,
        extract_bin_size=1000,
    )

    submit_stub.assert_called_once()
    kwargs = submit_stub.call_args.kwargs
    assert kwargs["display_name"] == "Nexus SOMA Randomized Extract - test_extract"
    assert kwargs["gcp_project"] == "test-project"
    assert kwargs["pipeline_kwargs"]["extract_configs"] == ["gs://bucket/extract_1.yaml"]
    assert kwargs["pipeline_kwargs"]["mark_finished_config"] == "gs://bucket/extract_1.yaml"
    assert url == "https://console.cloud.google.com/vertex-ai/pipelines/123"
