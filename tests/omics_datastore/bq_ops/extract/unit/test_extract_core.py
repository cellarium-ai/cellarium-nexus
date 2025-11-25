import typing

import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import extract as extract_module
from cellarium.nexus.shared import schemas as schemas_module


class FakeStorageReader:
    def __init__(self, rows: list[dict[str, typing.Any]]) -> None:
        self._rows = rows

    def rows(self, session: typing.Any) -> list[dict[str, typing.Any]]:  # noqa: D401
        return self._rows


class FakeStorageClient:
    def __init__(self) -> None:
        self.last_parent: str | None = None
        self.last_requested_session: typing.Any | None = None
        self.last_max_stream_count: int | None = None
        self._rows: list[dict[str, typing.Any]] = []
        self._session = type(
            "_Session", (), {"estimated_total_bytes_scanned": 0, "streams": [type("_S", (), {"name": "stream-1"})()]}
        )()

    def set_rows(self, rows: list[dict[str, typing.Any]]) -> None:
        self._rows = rows

    def create_read_session(self, *, parent: str, read_session: typing.Any, max_stream_count: int) -> typing.Any:
        self.last_parent = parent
        self.last_requested_session = read_session
        self.last_max_stream_count = max_stream_count
        return self._session

    def read_rows(self, name: str) -> FakeStorageReader:
        return FakeStorageReader(rows=self._rows)


def test_execute_query_delegates_to_client(bq_client: typing.Any) -> None:
    """
    Verify that execute_query delegates to the client's query and waits for result.

    :param bq_client: BigQuery client mock fixture
    """
    de = extract_module.DataExtractor(
        client=bq_client,
        project="p",
        dataset="d",
        extract_table_prefix="x_",
    )
    res = de.execute_query(sql="SELECT 1")
    assert bq_client.query_sql_recorder[-1] == "SELECT 1"
    assert res is None


def test_get_features_renders_and_builds_models(monkeypatch: pytest.MonkeyPatch, bq_client: typing.Any) -> None:
    """
    Validate that get_features renders the expected template and returns FeatureSchema models.

    :param monkeypatch: Pytest monkeypatch fixture
    :param bq_client: BigQuery client mock fixture
    """
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        return "-- get features"

    monkeypatch.setattr(extract_module.bq_sql, "render", _render)

    rows = [
        {"id": 1, "symbol": "A", "ensemble_id": "ens1"},
        {"id": 2, "symbol": "B", "ensemble_id": "ens2"},
    ]

    de = extract_module.DataExtractor(
        client=bq_client,
        project="proj",
        dataset="ds",
        extract_table_prefix="pref_",
    )
    monkeypatch.setattr(de, "execute_query", lambda sql: rows)

    feats = de.get_features()

    # ensure template used
    path, _td = render_calls[0]
    assert str(extract_module.GET_FEATURES_TEMPLATE) in path
    # ensure pydantic models built
    assert feats == [
        schemas_module.FeatureSchema(id=1, symbol="A", ensemble_id="ens1"),
        schemas_module.FeatureSchema(id=2, symbol="B", ensemble_id="ens2"),
    ]


def test_get_cells_in_bin_range_renders_with_params(monkeypatch: pytest.MonkeyPatch, bq_client: typing.Any) -> None:
    """
    Ensure get_cells_in_bin_range renders with start/end bins and optional select columns.

    :param monkeypatch: Pytest monkeypatch fixture
    :param bq_client: BigQuery client mock fixture
    """
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        return "-- get cells"

    monkeypatch.setattr(extract_module.bq_sql, "render", _render)

    de = extract_module.DataExtractor(
        client=bq_client,
        project="P",
        dataset="D",
        extract_table_prefix="X_",
    )
    # Return two rows
    monkeypatch.setattr(de, "execute_query", lambda sql: [{"id": 1}, {"id": 2}])

    res = de.get_cells_in_bin_range(start_bin=10, end_bin=11, obs_columns=["c1", "c2"])
    assert len(res) == 2

    _path, td = render_calls[0]
    assert str(extract_module.GET_CELLS_IN_BIN_RANGE_TEMPLATE) in _path
    assert td.data["project"] == "P"
    assert td.data["dataset"] == "D"
    assert td.data["extract_table_prefix"] == "X_"
    assert td.data["start_bin"] == 10
    assert td.data["end_bin"] == 11
    assert td.data["select_columns"] == ["c1", "c2"]


def test_get_matrix_data_builds_read_session(monkeypatch: pytest.MonkeyPatch, bq_client: typing.Any) -> None:
    """
    Confirm get_matrix_data builds a BigQuery Storage read session with proper fields and filter.

    :param monkeypatch: Pytest monkeypatch fixture
    :param bq_client: BigQuery client mock fixture
    """
    de = extract_module.DataExtractor(
        client=bq_client,
        project="proj",
        dataset="ds",
        extract_table_prefix="pref_",
    )

    fake_storage = FakeStorageClient()
    fake_storage.set_rows([{"cell_id": 100, "feature_data": [{"feature_id": 1, "raw_counts": 2.0}]}])
    # Inject fake BigQuery Storage client
    de.read_client = fake_storage  # type: ignore[assignment]

    out = de.get_matrix_data(start_bin=5, end_bin=6)

    # assertions about requested session
    assert fake_storage.last_parent == "projects/proj"
    assert fake_storage.last_max_stream_count == 1
    # verify selected fields & row restriction
    rs = fake_storage.last_requested_session
    assert rs.read_options.selected_fields == ["cell_id", "feature_data"]
    assert "extract_bin BETWEEN 5 AND 6" in rs.read_options.row_restriction

    # ensure returned rows match fake
    assert list(out) == [{"cell_id": 100, "feature_data": [{"feature_id": 1, "raw_counts": 2.0}]}]
