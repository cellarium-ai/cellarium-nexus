"""
Unit tests for SOMA workflows utilities.
"""

from unittest import mock

import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.cell_management.utils import exceptions, soma_workflows_utils
from cellarium.nexus.shared.schemas import component_configs


@pytest.fixture
def soma_dataset():
    """Create a mock SOMA dataset."""
    dataset = mock.Mock(spec=models.OmicsDataset)
    dataset.name = "test_soma_dataset"
    dataset.uri = "gs://bucket/soma_experiment"
    dataset.backend = models.OmicsDatasetBackend.TILEDB_SOMA
    return dataset


class TestGetTotalCellsInSoma:
    """Tests for get_total_cells_in_soma function."""

    def test_counts_cells_with_filters(self, soma_dataset):
        """Test that count_cells is called with correct filters."""
        filters = {"cell_type__in": ["T cell", "B cell"]}

        with mock.patch(
            "cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.TileDBSOMADataOperator"
        ) as mock_operator_cls:
            mock_operator = mock.Mock()
            mock_operator.count_cells.return_value = 1000
            mock_operator_cls.return_value = mock_operator

            result = soma_workflows_utils.get_total_cells_in_soma(omics_dataset=soma_dataset, filters=filters)

            mock_operator_cls.assert_called_once_with(experiment_uri="gs://bucket/soma_experiment")
            mock_operator.count_cells.assert_called_once_with(filter_statements=filters)
            assert result == 1000

    def test_raises_error_when_no_uri(self, soma_dataset):
        """Test that ValueError is raised when dataset has no URI."""
        soma_dataset.uri = None

        with pytest.raises(ValueError, match="has no URI configured"):
            soma_workflows_utils.get_total_cells_in_soma(omics_dataset=soma_dataset)


class TestComposeSomaExtractConfigs:
    """Tests for compose_soma_extract_configs function."""

    def test_invalid_range_size_raises_error(self, soma_dataset):
        """Test that ValueError is raised for invalid range_size."""
        with pytest.raises(ValueError, match="Range size must be greater than 0"):
            soma_workflows_utils.compose_soma_extract_configs(
                name="test_extract",
                creator_id=1,
                omics_dataset=soma_dataset,
                range_size=0,
                output_chunk_size=1000,
            )

    def test_no_uri_raises_error(self, soma_dataset):
        """Test that ValueError is raised when dataset has no URI."""
        soma_dataset.uri = None

        with pytest.raises(ValueError, match="has no URI configured"):
            soma_workflows_utils.compose_soma_extract_configs(
                name="test_extract",
                creator_id=1,
                omics_dataset=soma_dataset,
                range_size=1000,
                output_chunk_size=1000,
            )

    def test_zero_cells_raises_error(self, soma_dataset):
        """Test that ZeroCellsReturnedError is raised when no cells match."""
        with mock.patch(
            "cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.get_total_cells_in_soma",
            return_value=0,
        ):
            with pytest.raises(exceptions.ZeroCellsReturnedError):
                soma_workflows_utils.compose_soma_extract_configs(
                    name="test_extract",
                    creator_id=1,
                    omics_dataset=soma_dataset,
                    range_size=1000,
                    output_chunk_size=1000,
                )

    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.settings")
    def test_happy_path_creates_configs(self, mock_settings, soma_dataset):
        """Test that configs are created correctly."""
        mock_settings.BACKEND_PIPELINE_DIR = "backend"
        mock_settings.PIPELINE_DATA_EXTRACTS_DIR = "extracts"
        mock_settings.SITE_URL = "https://nexus.example.com"
        mock_settings.BUCKET_NAME_PRIVATE = "private-bucket"

        with mock.patch(
            "cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.get_total_cells_in_soma",
            return_value=5000,
        ):
            prepare_config, extract_configs = soma_workflows_utils.compose_soma_extract_configs(
                name="test_extract",
                creator_id=42,
                omics_dataset=soma_dataset,
                range_size=1000,
                output_chunk_size=1000,
                filters={"cell_type__eq": "T cell"},
                obs_columns=["cell_type", "tissue"],
            )

        # Verify prepare config
        assert isinstance(prepare_config, component_configs.SomaOpsPrepareExtract)
        assert prepare_config.extract_name == "test_extract"
        assert prepare_config.creator_id == 42
        assert prepare_config.experiment_uri == "gs://bucket/soma_experiment"
        assert prepare_config.range_size == 1000
        assert prepare_config.output_chunk_size == 1000
        assert prepare_config.filters == {"cell_type__eq": "T cell"}
        assert prepare_config.plan_path == "backend/extracts/test_extract/soma_extract_plan.json"

        # Verify extract configs (5000 cells / 1000 range_size = 5 ranges, all fit in one worker batch)
        assert len(extract_configs) == 1
        assert extract_configs[0].extract_name == "test_extract"
        assert extract_configs[0].range_indices == [0, 1, 2, 3, 4]
        assert extract_configs[0].obs_columns == ["cell_type", "tissue"]

    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.settings")
    def test_multiple_worker_batches(self, mock_settings, soma_dataset):
        """Test that multiple worker batches are created for large extracts."""
        mock_settings.BACKEND_PIPELINE_DIR = "backend"
        mock_settings.PIPELINE_DATA_EXTRACTS_DIR = "extracts"
        mock_settings.SITE_URL = "https://nexus.example.com"
        mock_settings.BUCKET_NAME_PRIVATE = "private-bucket"

        # 100 ranges with BINS_PER_WORKER=32 should create 4 worker batches
        with mock.patch(
            "cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.get_total_cells_in_soma",
            return_value=100_000,
        ):
            prepare_config, extract_configs = soma_workflows_utils.compose_soma_extract_configs(
                name="test_extract",
                creator_id=1,
                omics_dataset=soma_dataset,
                range_size=1000,  # 100 ranges
                output_chunk_size=1000,
            )

        # 100 ranges / 32 per worker = 4 worker batches (32, 32, 32, 4)
        assert len(extract_configs) == 4
        assert extract_configs[0].range_indices == list(range(0, 32))
        assert extract_configs[1].range_indices == list(range(32, 64))
        assert extract_configs[2].range_indices == list(range(64, 96))
        assert extract_configs[3].range_indices == list(range(96, 100))


class TestComposeAndDumpSomaConfigs:
    """Tests for compose_and_dump_soma_configs function."""

    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.workflows_configs")
    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.compose_soma_extract_configs")
    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.settings")
    def test_delegates_and_returns_paths(self, mock_settings, mock_compose, mock_workflows_configs, soma_dataset):
        """Test that configs are composed and dumped to GCS."""
        mock_settings.BUCKET_NAME_PRIVATE = "private-bucket"
        mock_settings.PIPELINE_CONFIGS_DIR = "configs"

        mock_prepare_config = mock.Mock()
        mock_extract_configs = [mock.Mock(), mock.Mock()]
        mock_compose.return_value = (mock_prepare_config, mock_extract_configs)

        mock_workflows_configs.dump_configs_to_bucket.side_effect = [
            ["gs://bucket/prepare_config.yaml"],
            ["gs://bucket/extract_1.yaml", "gs://bucket/extract_2.yaml"],
        ]

        prepare_path, extract_paths = soma_workflows_utils.compose_and_dump_soma_configs(
            name="test_extract",
            creator_id=1,
            omics_dataset=soma_dataset,
            range_size=1000,
            output_chunk_size=1000,
        )

        assert prepare_path == "gs://bucket/prepare_config.yaml"
        assert extract_paths == ["gs://bucket/extract_1.yaml", "gs://bucket/extract_2.yaml"]


class TestSubmitSomaExtractPipeline:
    """Tests for submit_soma_extract_pipeline function."""

    @mock.patch("cellarium.nexus.workflows.kubeflow.utils.job.submit_pipeline")
    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.compose_and_dump_soma_configs")
    @mock.patch("cellarium.nexus.backend.cell_management.utils.soma_workflows_utils.settings")
    def test_calls_submit_pipeline(self, mock_settings, mock_compose_and_dump, mock_submit, soma_dataset):
        """Test that pipeline is submitted with correct parameters."""
        mock_settings.GCP_PROJECT_ID = "test-project"
        mock_settings.PIPELINE_SERVICE_ACCOUNT = "sa@test.iam.gserviceaccount.com"
        mock_settings.PIPELINE_ROOT_PATH = "gs://bucket/pipeline-root"
        mock_settings.GCP_APPLICATION_BILLING_LABEL = "nexus"

        mock_compose_and_dump.return_value = (
            "gs://bucket/prepare_config.yaml",
            ["gs://bucket/extract_1.yaml"],
        )

        mock_submit.return_value = "https://console.cloud.google.com/vertex-ai/pipelines/123"

        result = soma_workflows_utils.submit_soma_extract_pipeline(
            name="test_extract",
            creator_id=1,
            omics_dataset=soma_dataset,
            range_size=1000,
            output_chunk_size=1000,
        )

        mock_submit.assert_called_once()
        call_kwargs = mock_submit.call_args[1]
        assert call_kwargs["display_name"] == "Nexus SOMA Extract - test_extract"
        assert call_kwargs["gcp_project"] == "test-project"
        assert call_kwargs["pipeline_kwargs"]["prepare_extract_config"] == "gs://bucket/prepare_config.yaml"
        assert call_kwargs["pipeline_kwargs"]["extract_configs"] == ["gs://bucket/extract_1.yaml"]
        assert result == "https://console.cloud.google.com/vertex-ai/pipelines/123"
