from __future__ import annotations

from pathlib import Path
from typing import Callable

import pytest

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.backend.curriculum import models as curriculum_models


@pytest.fixture()
def default_dataset() -> models.BigQueryDataset:
    """
    Create default BigQuery dataset row for backend tests.

    :return: BigQueryDataset instance
    """
    return models.BigQueryDataset.objects.create(name="default", description="Default dataset")


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
