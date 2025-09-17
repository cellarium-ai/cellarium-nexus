import concurrent.futures as cf
from typing import Iterator

import numpy as np
import pytest

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
    Patch ProcessPoolExecutor in the target module to use ThreadPoolExecutor to avoid multiprocessing in tests.

    :param monkeypatch: Pytest monkeypatch fixture

    :return: None
    """
    monkeypatch.setattr(create_ingest_files.concurrency, "ProcessPoolExecutor", cf.ThreadPoolExecutor)
    yield
