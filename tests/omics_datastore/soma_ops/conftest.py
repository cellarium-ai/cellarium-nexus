"""
Shared fixtures and mock classes for soma_ops tests.

This module provides reusable mock classes for testing SOMA and AnnData operations.
"""

from pathlib import Path
from typing import Callable
from unittest.mock import MagicMock

import pandas as pd
import scipy.sparse as sp


class FakeAnnData:
    """Mock AnnData object for testing."""

    def __init__(self, X: object = None, obs: object = None, var: object = None, n_obs: int = 0) -> None:
        self.X = X
        self.obs = obs
        self.var = var
        self.n_obs = n_obs
        self.file = MagicMock()

    def __getitem__(self, indices: object) -> "FakeAnnData":
        if isinstance(indices, list):
            return FakeAnnData(n_obs=len(indices))
        return FakeAnnData(n_obs=1)

    def copy(self) -> "FakeAnnData":
        return FakeAnnData(n_obs=self.n_obs)

    def write_h5ad(self, path: Path, compression: str | None = None) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    def write_zarr(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.mkdir(exist_ok=True)


class FakeFuture:
    """Mock concurrent.futures.Future for testing."""

    def __init__(self, result: object = None, exception: Exception | None = None) -> None:
        self._result = result
        self._exception = exception

    def result(self) -> object:
        if self._exception:
            raise self._exception
        return self._result


class FakeExecutor:
    """Mock ProcessPoolExecutor for testing."""

    def __init__(self, max_workers: int, worker_fn: Callable | None = None) -> None:
        self.max_workers = max_workers
        self.worker_fn = worker_fn
        self.futures: list[FakeFuture] = []
        self.calls: list[tuple[object, tuple, dict]] = []

    def __enter__(self) -> "FakeExecutor":
        return self

    def __exit__(self, *args: object) -> None:
        pass

    def submit(self, fn: object, *args: object, **kwargs: object) -> FakeFuture:
        self.calls.append((fn, args, kwargs))
        if self.worker_fn:
            result = self.worker_fn(*args, **kwargs)
            future = FakeFuture(result=result)
        else:
            # Execute synchronously by default
            result = fn(*args, **kwargs)  # type: ignore[operator]
            future = FakeFuture(result=result)
        self.futures.append(future)
        return future


class FakeSomaPandasTable:
    """Mock SOMA table that returns a pandas DataFrame."""

    def __init__(self, data: pd.DataFrame) -> None:
        self.data = data

    def to_pandas(self) -> pd.DataFrame:
        return self.data.copy()


class FakeSomaQuery:
    """Mock SOMA query that returns a table."""

    def __init__(self, data: pd.DataFrame) -> None:
        self.data = data

    def concat(self) -> FakeSomaPandasTable:
        return FakeSomaPandasTable(self.data)


class FakeSomaObs:
    """Mock SOMA obs dataframe."""

    def __init__(self, data: pd.DataFrame, expected_filter: str | None = None) -> None:
        self.data = data
        self.expected_filter = expected_filter

    def read(
        self,
        coords: tuple[object] | None = None,
        *,
        column_names: list[str] | None = None,
        value_filter: str | None = None
    ) -> FakeSomaQuery:
        # Validate column_names if provided
        if column_names is not None:
            assert column_names == ["soma_joinid"]
        # Validate expected_filter if set
        if self.expected_filter is not None:
            assert value_filter == self.expected_filter or value_filter is None
        return FakeSomaQuery(self.data)


class FakeSomaXTable:
    """Mock SOMA X sparse matrix table."""

    def __init__(self, matrix: sp.coo_matrix) -> None:
        self.matrix = matrix

    def to_scipy(self, compress: bool = False) -> sp.coo_matrix:
        return self.matrix


class FakeSomaXQuery:
    """Mock SOMA X query."""

    def __init__(self, matrix: sp.coo_matrix) -> None:
        self.matrix = matrix

    def concat(self) -> FakeSomaXTable:
        return FakeSomaXTable(self.matrix)


class FakeSomaXData:
    """Mock SOMA X data layer."""

    def __init__(self, matrix: sp.coo_matrix) -> None:
        self.matrix = matrix

    def read(self, coords: tuple[object, object]) -> FakeSomaXQuery:
        return FakeSomaXQuery(self.matrix)


class FakeSomaVar:
    """Mock SOMA var dataframe."""

    def __init__(self, data: pd.DataFrame) -> None:
        self.data = data

    def read(self) -> FakeSomaQuery:
        return FakeSomaQuery(self.data)


class FakeSomaMeasurement:
    """Mock SOMA measurement."""

    def __init__(self, var_data: pd.DataFrame, x_matrix: sp.coo_matrix, x_layer: str = "X") -> None:
        self.var = FakeSomaVar(var_data)
        self.X = {x_layer: FakeSomaXData(x_matrix)}


class FakeSomaExperiment:
    """
    Mock SOMA experiment.

    Supports both simple (obs-only) and complex (obs + var + X) configurations.
    """

    def __init__(
        self,
        obs_data: pd.DataFrame,
        var_data: pd.DataFrame | None = None,
        x_matrix: sp.coo_matrix | None = None,
        measurement_name: str = "RNA",
        x_layer: str = "X",
        expected_filter: str | None = None,
    ) -> None:
        self.obs = FakeSomaObs(obs_data, expected_filter)
        if var_data is not None and x_matrix is not None:
            self.ms = {measurement_name: FakeSomaMeasurement(var_data, x_matrix, x_layer)}
        else:
            self.ms = {}

    def __enter__(self) -> "FakeSomaExperiment":
        return self

    def __exit__(self, *args: object) -> None:
        pass
