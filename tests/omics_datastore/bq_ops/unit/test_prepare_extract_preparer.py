import typing

import pytest
from google.cloud import bigquery

from cellarium.nexus.omics_datastore.bq_ops import constants as extract_constants
from cellarium.nexus.omics_datastore.bq_ops.extract import prepare_extract as prepare_extract_module
from cellarium.nexus.shared import schemas as schemas_module


class TableObj:
    def __init__(self) -> None:
        self.expires = None


class ClientWithExpiration:
    def __init__(self) -> None:
        self.get_table_calls: list[str] = []
        self.update_table_calls: list[tuple[TableObj, list[str]]] = []

    def get_table(self, table_id: str) -> TableObj:
        self.get_table_calls.append(table_id)
        return TableObj()

    def update_table(self, table: TableObj, fields: list[str]) -> None:
        self.update_table_calls.append((table, fields))


class ClientForFeatureLoad:
    def __init__(self) -> None:
        self.loaded: list[tuple[str, bigquery.LoadJobConfig]] = []
        self._table = TableObj()
        self.update_table_calls: list[list[str]] = []

    def load_table_from_file(self, file_obj: typing.Any, table_id: str, job_config: bigquery.LoadJobConfig):
        # Record table id and job_config; return object with result() and output_rows
        self.loaded.append((table_id, job_config))

        class _Job:
            def __init__(self) -> None:
                self.output_rows = 5

            def result(self) -> None:
                return None

        return _Job()

    def get_table(self, table_id: str) -> TableObj:
        return self._table

    def update_table(self, table: TableObj, fields: list[str]) -> None:
        self.update_table_calls.append(fields)


class ExecPreparer(prepare_extract_module.ExtractTablePreparer):
    """
    Preparer subclass that records executed SQL statements for assertions.
    """

    def __init__(self, *, client: typing.Any, project: str, dataset: str, extract_table_prefix: str) -> None:
        super().__init__(client=client, project=project, dataset=dataset, extract_table_prefix=extract_table_prefix)
        self.executed_sql: list[str] = []

    def execute_query(self, sql: str) -> typing.Any:  # type: ignore[override]
        self.executed_sql.append(sql)

        class _Job:
            def result(self) -> None:
                return None

        return _Job()


def test_execute_query_calls_client_query_and_result(bq_client: typing.Any) -> None:
    """
    Verify that execute_query delegates to client's query and waits for result.

    :param bq_client: BigQuery client mock fixture
    """
    preparer = prepare_extract_module.ExtractTablePreparer(
        client=bq_client,
        project="p",
        dataset="d",
        extract_table_prefix="x_",
    )

    res = preparer.execute_query("SELECT 1")
    # BQClientMock records SQL; result() is available on job mock
    assert bq_client.query_sql_recorder[-1] == "SELECT 1"
    # result() returns None in our mock; execute_query returns that value
    assert res is None


def test_set_table_expiration_sets_expires_and_updates_table(freeze_time: typing.Callable[[typing.Any], None]) -> None:
    """
    Ensure set_table_expiration sets UTC expiry and updates the table.

    :param freeze_time: Fixture to freeze the current time
    """
    import datetime as dt

    client = ClientWithExpiration()
    preparer = prepare_extract_module.ExtractTablePreparer(
        client=client, project="p", dataset="d", extract_table_prefix="x_"
    )

    # Freeze to a known instant with tz awareness
    freeze_time(dt.datetime(2025, 1, 1, 0, 0, 0, tzinfo=dt.UTC))

    table_id = "p.d.t"
    preparer.set_table_expiration(table_id=table_id, expiration_delta=dt.timedelta(hours=3))

    assert client.get_table_calls == [table_id]
    assert len(client.update_table_calls) == 1
    (table_obj, fields) = client.update_table_calls[0]
    assert fields == ["expires"]
    assert table_obj.expires is not None
    assert table_obj.expires.tzinfo is not None


def test_prepare_feature_table_loads_csv_and_sets_expiration(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Validate that prepare_feature_table loads CSV to BQ and sets expiration.

    :param monkeypatch: Pytest monkeypatch fixture
    """
    # Build some feature schemas
    feats = [
        schemas_module.FeatureSchema(id=1, symbol="A", ensemble_id="ens1"),
        schemas_module.FeatureSchema(id=2, symbol="B", ensemble_id="ens2"),
    ]

    client = ClientForFeatureLoad()

    preparer = prepare_extract_module.ExtractTablePreparer(
        client=client, project="proj", dataset="ds", extract_table_prefix="pref_"
    )

    preparer.prepare_feature_table(feats)

    # One load issued
    assert len(client.loaded) == 1
    (table_id, job_config) = client.loaded[0]
    assert table_id == f"proj.ds.pref_{extract_constants.BQ_EXTRACT_FEATURE_INFO_TABLE_NAME}"
    assert job_config.source_format == bigquery.SourceFormat.CSV
    assert job_config.write_disposition == bigquery.WriteDisposition.WRITE_TRUNCATE

    # Expiration update called with fields=["expires"]
    assert client.update_table_calls and client.update_table_calls[-1] == ["expires"]


def test_prepare_cell_info_renders_and_executes_in_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Check that prepare_cell_info renders and executes templates in order.

    The expected order is randomized → binned → drop. Also assert TemplateData
    fields and that the final cell info table receives a 14-day expiration.

    :param monkeypatch: Pytest monkeypatch fixture
    """
    # Capture render calls and template data
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        # return SQL including path to allow order assertions
        return f"-- {path}\nSELECT 1;"

    monkeypatch.setattr(prepare_extract_module.bq_sql, "render", _render)

    exp_calls: list[tuple[str, typing.Any]] = []

    def _set_expiration(*, table_id: str, expiration_delta: typing.Any) -> None:
        exp_calls.append((table_id, expiration_delta))

    prep = ExecPreparer(client=object(), project="p", dataset="d", extract_table_prefix="x_")
    exp_calls: list[tuple[str, typing.Any]] = []

    def _set_expiration(*, table_id: str, expiration_delta: typing.Any) -> None:
        exp_calls.append((table_id, expiration_delta))

    monkeypatch.setattr(
        prep,
        "set_table_expiration",
        lambda table_id, expiration_delta: _set_expiration(table_id=table_id, expiration_delta=expiration_delta),
    )
    monkeypatch.setattr(
        prep,
        "set_table_expiration",
        lambda table_id, expiration_delta: _set_expiration(table_id=table_id, expiration_delta=expiration_delta),
    )

    prep.prepare_cell_info(
        extract_bin_size=100,
        random_seed_offset=7,
        partition_bin_count=123,
        partition_size=9,
        extract_bin_keys=["a", "b"],
        filters={"k__eq": "v"},
        obs_columns=["c1", "c2"],
        metadata_extra_columns=["m1"],
    )

    # Three executes: randomized, binned, drop
    assert len(prep.executed_sql) == 3
    # Ensure the renders were invoked with expected templates
    template_paths = [p for (p, _d) in render_calls]
    assert str(prepare_extract_module.CELL_INFO_RAND_TEMPLATE) in template_paths[0]
    assert str(prepare_extract_module.CELL_INFO_TEMPLATE) in template_paths[1]
    assert str(prepare_extract_module.DROP_CELL_INFO_RAND_TEMPLATE) in template_paths[2]

    # Expiration set for cell info table with 14 days
    assert len(exp_calls) == 1
    table_id, delta = exp_calls[0]
    assert table_id == f"p.d.x_{extract_constants.BQ_EXTRACT_CELL_INFO_TABLE_NAME}"
    import datetime as dt

    assert isinstance(delta, dt.timedelta) and delta.days == 14

    # Verify TemplateData fields populated as expected for the binned table call (second render)
    _, td = render_calls[1]
    assert td.data["random_seed_offset"] == 7
    assert td.data["partition_bin_count"] == 123
    assert td.data["partition_size"] == 9
    assert td.data["extract_bin_size"] == 100
    assert td.data["extract_bin_keys"] == ["a", "b"]
    assert td.data["metadata_columns"] == ["m1"]
    assert td.data["select_columns"] == ["c1", "c2"]
    assert td.data["filter_statements"] == {"k__eq": "v"}


def test_prepare_count_matrix_renders_and_sets_expiration(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    Verify that prepare_count_matrix renders SQL, executes it, and sets expiry.

    :param monkeypatch: Pytest monkeypatch fixture
    """
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        return f"-- {path}\nSELECT 2;"

    monkeypatch.setattr(prepare_extract_module.bq_sql, "render", _render)

    exp_calls: list[tuple[str, typing.Any]] = []

    def _set_expiration(*, table_id: str, expiration_delta: typing.Any) -> None:
        exp_calls.append((table_id, expiration_delta))

    prep = ExecPreparer(client=object(), project="p", dataset="d", extract_table_prefix="x_")
    monkeypatch.setattr(
        prep,
        "set_table_expiration",
        lambda table_id, expiration_delta: _set_expiration(table_id=table_id, expiration_delta=expiration_delta),
    )

    prep.prepare_count_matrix(partition_bin_count=456, partition_size=12)

    # One execute and one expiration call
    assert len(prep.executed_sql) == 1
    assert len(exp_calls) == 1
    table_id, delta = exp_calls[0]
    assert table_id == f"p.d.x_{extract_constants.BQ_EXTRACT_MATRIX_COO_TABLE_NAME}"

    # Verify TemplateData fields
    _path, td = render_calls[0]
    assert td.data["partition_bin_count"] == 456
    assert td.data["partition_size"] == 12


@pytest.mark.parametrize(
    "extract_bin_keys,metadata_extra_columns,filters,obs_columns,expect_select",
    [
        (None, None, {"a__eq": 1}, None, False),
        (["k1"], None, {"a__in": [1, 2]}, None, False),
        (None, ["m1"], {"a__eq": 1}, ["c1"], True),
        (["k1", "k2"], ["m1", "m2"], {"a__in": [3, 4]}, ["c1", "c2"], True),
    ],
)
def test_prepare_cell_info_template_data_variants(
    monkeypatch: pytest.MonkeyPatch,
    extract_bin_keys: list[str] | None,
    metadata_extra_columns: list[str] | None,
    filters: dict[str, typing.Any],
    obs_columns: list[str] | None,
    expect_select: bool,
) -> None:
    """
    Validate TemplateData population across parameter variants for prepare_cell_info.

    :param monkeypatch: Pytest monkeypatch fixture
    :param extract_bin_keys: Keys for binning or None
    :param metadata_extra_columns: Extra metadata columns or None
    :param filters: Filter dictionary using '<col>__<op>' format
    :param obs_columns: Observation columns or None
    :param expect_select: Whether select_columns is expected to be present
    """
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        return f"-- {path}\nSELECT 1;"

    monkeypatch.setattr(prepare_extract_module.bq_sql, "render", _render)

    prep = ExecPreparer(client=object(), project="p", dataset="d", extract_table_prefix="x_")
    exp_calls_local: list[tuple[str, typing.Any]] = []

    def _set_exp_local(*, table_id: str, expiration_delta: typing.Any) -> None:
        exp_calls_local.append((table_id, expiration_delta))

    monkeypatch.setattr(
        prep,
        "set_table_expiration",
        lambda table_id, expiration_delta: _set_exp_local(table_id=table_id, expiration_delta=expiration_delta),
    )

    # Only need to run until binned step to observe full TemplateData; use the method directly
    prep.prepare_cell_info(
        extract_bin_size=None,
        random_seed_offset=0,
        partition_bin_count=10,
        partition_size=2,
        extract_bin_keys=extract_bin_keys,
        filters=filters,
        obs_columns=obs_columns,
        metadata_extra_columns=metadata_extra_columns,
    )

    # Second render corresponds to the binned table
    _, td = render_calls[1]
    # Common assertions
    assert td.data["partition_bin_count"] == 10
    assert td.data["partition_size"] == 2
    assert td.data["extract_bin_size"] is None
    assert td.data["extract_bin_keys"] == extract_bin_keys
    # metadata columns echoed as provided
    assert td.data["metadata_columns"] == metadata_extra_columns
    # filters stored under filter_statements
    assert td.data["filter_statements"] == filters
    # select presence based on obs_columns
    if expect_select:
        assert td.data["select_columns"] == obs_columns
    else:
        assert "select_columns" not in td.data


def test_prepare_count_matrix_defaults(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify prepare_count_matrix populates TemplateData defaults when no overrides are passed.

    :param monkeypatch: Pytest monkeypatch fixture
    """
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        return f"-- {path}\nSELECT 1;"

    monkeypatch.setattr(prepare_extract_module.bq_sql, "render", _render)

    prep = ExecPreparer(client=object(), project="p", dataset="d", extract_table_prefix="x_")
    exp_calls: list[tuple[str, typing.Any]] = []

    def _set_expiration(*, table_id: str, expiration_delta: typing.Any) -> None:
        exp_calls.append((table_id, expiration_delta))

    monkeypatch.setattr(
        prep,
        "set_table_expiration",
        lambda table_id, expiration_delta: _set_expiration(table_id=table_id, expiration_delta=expiration_delta),
    )
    prep.prepare_count_matrix()

    # One render call, verify defaults from function signature
    assert len(render_calls) == 1
    _path, td = render_calls[0]
    assert td.data["partition_bin_count"] == 40000
    assert td.data["partition_size"] == 10
