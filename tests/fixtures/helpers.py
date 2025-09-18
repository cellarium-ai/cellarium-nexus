import concurrent.futures as cf
from typing import Iterator

import numpy as np
import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import extract as extract_module
from cellarium.nexus.omics_datastore.bq_ops.ingest import create_ingest_files


@pytest.fixture(autouse=True)
def rng_seed() -> None:
    """
    Set deterministic random seed for all tests that load this plugin.

    :return: None
    """
    np.random.seed(0)


@pytest.fixture()
def patched_process_pool_executor(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """
    Patch ProcessPoolExecutor globally to use ThreadPoolExecutor to avoid multiprocessing in tests.

    :param monkeypatch: Pytest monkeypatch fixture

    :return: None
    """
    # Patch on the concurrent.futures module so any module aliasing it (e.g., "import concurrent.futures as concurrency")
    # sees the patched ProcessPoolExecutor.
    monkeypatch.setattr(cf, "ProcessPoolExecutor", cf.ThreadPoolExecutor)
    yield
