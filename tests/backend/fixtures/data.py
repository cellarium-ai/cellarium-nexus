from __future__ import annotations

from pathlib import Path
from typing import Callable

import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.curriculum import models as curriculum_models


@pytest.fixture()
def feature_schema_factory() -> Callable:
    """
    Provide factory to create FeatureSchema instances with features.

    :return: Callable that creates FeatureSchema with associated FeatureInfo records
    """

    def _create(name: str, features: list[tuple[str, str]] | None = None) -> models.FeatureSchema:
        """
        Create a FeatureSchema with optional features.

        :param name: Name for the schema
        :param features: Optional list of (ensemble_id, symbol) tuples

        :return: Created FeatureSchema instance
        """
        schema = models.FeatureSchema.objects.create(name=name)
        if features:
            for ensemble_id, symbol in features:
                models.FeatureInfo.objects.create(
                    ensemble_id=ensemble_id,
                    symbol=symbol,
                    feature_schema=schema,
                )
        return schema

    return _create


@pytest.fixture()
def default_dataset() -> models.OmicsDataset:
    """
    Create default omics dataset row for backend tests.

    :return: OmicsDataset instance
    """
    return models.OmicsDataset.objects.create(name="default", description="Default dataset")


@pytest.fixture()
def curriculum_factory(admin_user: curriculum_models.UserModel) -> Callable[[str], curriculum_models.Curriculum]:
    """
    Provide factory to create curriculum instances bound to admin user.

    :param admin_user: Superuser fixture for ownership

    :return: Callable that persists ``Curriculum`` records with defaults
    """

    def _create(name: str, **overrides: object) -> curriculum_models.Curriculum:
        payload: dict[str, object] = {
            "name": name,
            "creator": admin_user,
            "cell_count": overrides.pop("cell_count", 1000),
            "extract_bin_size": overrides.pop("extract_bin_size", 500),
            "extract_bin_count": overrides.pop("extract_bin_count", 2),
            "filters_json": overrides.pop("filters_json", {"organism": "Homo sapiens"}),
            "status": overrides.pop(
                "status",
                curriculum_models.Curriculum.STATUS_PREPARE,
            ),
        }
        payload.update(overrides)
        curriculum = curriculum_models.Curriculum.objects.create(**payload)
        return curriculum

    return _create


@pytest.fixture()
def sample_csv_file(tmp_path: Path) -> Path:
    """
    Create a valid CSV file with feature data for testing.

    :param tmp_path: Temporary directory path

    :return: Path to CSV file
    """
    csv_path = tmp_path / "features.csv"
    csv_path.write_text(
        "ensemble_id,symbol\n" "ENSG00000141510,TP53\n" "ENSG00000157764,BRAF\n" "ENSG00000146648,EGFR\n"
    )
    return csv_path


@pytest.fixture()
def invalid_csv_file(tmp_path: Path) -> Path:
    """
    Create an invalid CSV file (missing required columns).

    :param tmp_path: Temporary directory path

    :return: Path to invalid CSV file
    """
    csv_path = tmp_path / "invalid.csv"
    csv_path.write_text("wrong_column,another_column\n" "value1,value2\n")
    return csv_path


@pytest.fixture()
def ingest_config_paths(tmp_path: Path) -> dict[str, Path]:
    """
    Prepare dummy Kubeflow config files for ingest pipelines.

    :param tmp_path: Base temporary directory fixture

    :return: Mapping with ``create_configs`` and ``ingest_config`` path references
    """

    configs_dir = tmp_path / "configs"
    configs_dir.mkdir()

    create_path = configs_dir / "create_ingest_files.yaml"
    create_path.write_text("kind: CreateIngestFilesConfig\n")

    ingest_path = configs_dir / "ingest_config.yaml"
    ingest_path.write_text("kind: IngestFilesConfig\n")

    return {
        "create_configs": create_path,
        "ingest_config": ingest_path,
    }
