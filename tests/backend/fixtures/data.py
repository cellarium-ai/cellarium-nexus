import pytest

from cellarium.nexus.backend.cell_management import models


@pytest.fixture()
def default_dataset() -> models.BigQueryDataset:
    """
    Create default BigQuery dataset row for backend tests.

    :return: BigQueryDataset instance
    """
    return models.BigQueryDataset.objects.create(name="default", description="Default dataset")
