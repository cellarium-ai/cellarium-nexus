import typing

import pytest
from google.api_core import exceptions as gapi_exceptions
from google.cloud import bigquery

from cellarium.nexus.omics_datastore.bq_ops import constants as ingest_constants
from cellarium.nexus.omics_datastore.bq_ops import exceptions as ingest_exceptions
from cellarium.nexus.omics_datastore.bq_ops.ingest import ingest_data_to_bigquery


def _record_create_staging_table() -> tuple[list[dict[str, typing.Any]], typing.Callable[..., None]]:
    calls: list[dict[str, typing.Any]] = []

    def _recorder(
        *,
        client: bigquery.Client,
        project: str,
        dataset: str,
        base_table_name: str,
        staging_suffix: str,
        schema: list[bigquery.SchemaField],
        clustering_fields: list[str] | None = None,
    ) -> None:
        calls.append(
            {
                "project": project,
                "dataset": dataset,
                "base_table_name": base_table_name,
                "staging_suffix": staging_suffix,
                "schema": schema,
                "clustering_fields": clustering_fields,
            }
        )

    return calls, _recorder


@pytest.fixture()
def patched_bq_client(monkeypatch: pytest.MonkeyPatch, bq_client: typing.Any) -> typing.Any:
    """
    Patch BigQuery client initializer to return our fixture-provided mock.

    :param monkeypatch: Pytest monkeypatch fixture

    :raise: None

    :return: Mock BigQuery client used by tests
    """
    monkeypatch.setattr(
        ingest_data_to_bigquery,
        "initialize_bigquery_client",
        lambda project_id, dataset: bq_client,
    )
    return bq_client


def test_context_happy_path_commits_and_cleans_up(
    monkeypatch: pytest.MonkeyPatch,
    patched_bq_client: typing.Any,
    freeze_time: typing.Callable[[typing.Any], None],
    freeze_uuid: typing.Callable[[str], None],
) -> None:
    """
    Exercise the happy path: staging creation → data loads → context exit → commit → cleanup.

    :param monkeypatch: Pytest monkeypatch fixture
    :param patched_bq_client: Mock BigQuery client
    :param freeze_time: Frozen time setter fixture
    :param freeze_uuid: Frozen uuid setter fixture

    :raise: AssertionError

    :return: None
    """
    # Make suffix deterministic
    import datetime as dt

    freeze_time(dt.datetime(2025, 1, 1, 0, 0, 0, tzinfo=dt.UTC))
    freeze_uuid("deadbeefcafebabe0123456789abcd0")
    expected_suffix = "staging_20250101000000_deadbe"

    # Record staging table creation calls
    calls, recorder = _record_create_staging_table()
    monkeypatch.setattr(ingest_data_to_bigquery, "create_staging_table", recorder)

    # Increase output rows for a more realistic job state
    patched_bq_client.load_job_rows = 123

    with ingest_data_to_bigquery.bigquery_ingest_context(
        project_id="proj",
        dataset="ds",
        gcs_bucket_name="bucket",
        gcs_stage_dir="stage",
    ):
        # No-op work inside context
        pass

    # Verify staging table creation
    assert len(calls) == 4
    base_names = {
        ingest_constants.BQ_INGEST_TABLE_NAME,
        ingest_constants.BQ_CELL_INFO_TABLE_NAME,
        ingest_constants.BQ_FEATURE_INFO_TABLE_NAME,
        ingest_constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME,
    }
    seen_bases = {c["base_table_name"] for c in calls}
    assert seen_bases == base_names
    # Clustering only for raw counts
    for c in calls:
        if c["base_table_name"] == ingest_constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME:
            assert c["clustering_fields"] == ["cell_id"]
        else:
            assert c["clustering_fields"] is None
        assert c["staging_suffix"] == expected_suffix

    # Verify loads issued to BQ
    # Expect four load calls; we only have the last one recorded by fixture,
    # so assert that at least the last call matches expected format.
    assert patched_bq_client.last_load_uri is not None
    assert patched_bq_client.last_load_table_id is not None
    assert patched_bq_client.last_load_job_config is not None
    assert patched_bq_client.last_load_uri.startswith("gs://bucket/stage/")
    assert patched_bq_client.last_load_table_id.startswith("proj.ds.")

    # Verify commit issued (multi-statement SQL)
    assert len(patched_bq_client.query_sql_recorder) == 1
    sql = patched_bq_client.query_sql_recorder[0]
    assert "BEGIN TRANSACTION;" in sql
    assert "COMMIT TRANSACTION;" in sql
    # Ensure inserts for each table are present with deterministic staging names
    for base in base_names:
        assert f"INSERT INTO `ds.{base}`" in sql
        assert f"FROM `ds.{base}_{expected_suffix}`" in sql

    # Verify cleanup (all staging tables deleted)
    # 4 tables expected
    assert len(patched_bq_client.delete_table_recorder) == 4
    for base in base_names:
        assert f"proj.ds.{base}_{expected_suffix}" in patched_bq_client.delete_table_recorder


def test_context_load_failure_raises_and_skips_commit(
    monkeypatch: pytest.MonkeyPatch,
    patched_bq_client: typing.Any,
    freeze_time: typing.Callable[[typing.Any], None],
    freeze_uuid: typing.Callable[[str], None],
) -> None:
    """
    Expect BigQueryLoadError when a staging load fails; commit should not run and cleanup should still happen.

    :param monkeypatch: Pytest monkeypatch fixture
    :param patched_bq_client: Mock BigQuery client
    :param freeze_time: Frozen time setter fixture
    :param freeze_uuid: Frozen uuid setter fixture

    :raise: AssertionError

    :return: None
    """
    import datetime as dt

    freeze_time(dt.datetime(2025, 1, 1, 0, 0, 0, tzinfo=dt.UTC))
    freeze_uuid("deadbeefcafebabe0123456789abcd0")

    # No-op staging creation
    _calls, recorder = _record_create_staging_table()
    monkeypatch.setattr(ingest_data_to_bigquery, "create_staging_table", recorder)

    # Fail the first load call without engaging tenacity retries
    def _failing_load(**kwargs: typing.Any) -> None:  # type: ignore[no-untyped-def]
        raise RuntimeError("boom")

    monkeypatch.setattr(ingest_data_to_bigquery, "perform_load_table_from_gcs", _failing_load)

    with pytest.raises(ingest_exceptions.BigQueryLoadError):
        with ingest_data_to_bigquery.bigquery_ingest_context(
            project_id="proj",
            dataset="ds",
            gcs_bucket_name="bucket",
            gcs_stage_dir="stage",
        ):
            pass

    # Commit skipped
    assert patched_bq_client.query_sql_recorder == []
    # Cleanup still executed (4 attempts)
    assert len(patched_bq_client.delete_table_recorder) == 4


def test_context_commit_failure_raises_and_cleans_up(
    monkeypatch: pytest.MonkeyPatch,
    patched_bq_client: typing.Any,
    freeze_time: typing.Callable[[typing.Any], None],
    freeze_uuid: typing.Callable[[str], None],
) -> None:
    """
    Expect BigQueryCommitError when commit step fails; staging cleanup still occurs.

    :param monkeypatch: Pytest monkeypatch fixture
    :param patched_bq_client: Mock BigQuery client
    :param freeze_time: Frozen time setter fixture
    :param freeze_uuid: Frozen uuid setter fixture

    :raise: AssertionError

    :return: None
    """
    import datetime as dt

    freeze_time(dt.datetime(2025, 1, 1, 0, 0, 0, tzinfo=dt.UTC))
    freeze_uuid("deadbeefcafebabe0123456789abcd0")

    # No-op staging creation
    _calls, recorder = _record_create_staging_table()
    monkeypatch.setattr(ingest_data_to_bigquery, "create_staging_table", recorder)

    # Make commit fail by raising from query
    def _failing_query(sql: str) -> typing.Any:
        raise RuntimeError("commit failed")

    monkeypatch.setattr(patched_bq_client, "query", _failing_query)

    with pytest.raises(ingest_exceptions.BigQueryCommitError):
        with ingest_data_to_bigquery.bigquery_ingest_context(
            project_id="proj",
            dataset="ds",
            gcs_bucket_name="bucket",
            gcs_stage_dir="stage",
        ):
            pass

    # Cleanup still executed (4 attempts)
    assert len(patched_bq_client.delete_table_recorder) == 4


def test_perform_load_table_from_gcs_retries_then_succeeds(
    monkeypatch: pytest.MonkeyPatch,
    patched_bq_client: typing.Any,
) -> None:
    """
    Validate retry behavior for perform_load_table_from_gcs.

    :param monkeypatch: Pytest monkeypatch fixture
    :param patched_bq_client: Mock BigQuery client

    :raise: AssertionError

    :return: None
    """
    # Make tenacity's backoff instant
    monkeypatch.setattr("time.sleep", lambda s: None)

    counter = {"n": 0}

    def _flaky_load(**kwargs: typing.Any) -> None:
        counter["n"] += 1
        if counter["n"] < 3:
            raise RuntimeError("transient")
        return None

    monkeypatch.setattr(ingest_data_to_bigquery, "load_table_from_gcs", _flaky_load)

    ingest_data_to_bigquery.perform_load_table_from_gcs(
        client=patched_bq_client,
        project="proj",
        dataset="ds",
        table_name="table",
        file_pattern="*.avro",
        file_format=bigquery.SourceFormat.AVRO,
        gcs_bucket_name="bucket",
        gcs_stage_dir="stage",
    )

    assert counter["n"] == 3


def test_get_schema_for_table_and_unknown_raises() -> None:
    """
    Check schema mapping for known tables and error for unknown tables.

    :raise: AssertionError

    :return: None
    """
    for base in (
        ingest_constants.BQ_INGEST_TABLE_NAME,
        ingest_constants.BQ_CELL_INFO_TABLE_NAME,
        ingest_constants.BQ_FEATURE_INFO_TABLE_NAME,
        ingest_constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME,
    ):
        schema = ingest_data_to_bigquery.get_schema_for_table(base)
        assert isinstance(schema, list)
        assert len(schema) > 0
        assert all(isinstance(f, bigquery.SchemaField) for f in schema)

    with pytest.raises(ingest_exceptions.BigQuerySchemaError):
        _ = ingest_data_to_bigquery.get_schema_for_table("unknown_table")


def test_load_table_from_gcs_builds_uri_and_table_id(patched_bq_client: typing.Any) -> None:
    """
    Verify URI and table_id composition for load_table_from_gcs.

    :param patched_bq_client: Mock BigQuery client

    :raise: AssertionError

    :return: None
    """
    ingest_data_to_bigquery.load_table_from_gcs(
        client=patched_bq_client,
        project="proj",
        dataset="ds",
        table_name="t1",
        file_pattern="file-*.avro",
        file_format=bigquery.SourceFormat.AVRO,
        gcs_bucket_name="b",
        gcs_stage_dir="s",
    )

    assert patched_bq_client.last_load_uri == "gs://b/s/file-*.avro"
    assert patched_bq_client.last_load_table_id == "proj.ds.t1"
    assert patched_bq_client.last_load_job_config is not None
    assert patched_bq_client.last_load_job_config.source_format == bigquery.SourceFormat.AVRO


def test_cleanup_staging_tables_handles_notfound(
    monkeypatch: pytest.MonkeyPatch,
    patched_bq_client: typing.Any,
) -> None:
    """
    Ensure cleanup_staging_tables handles NotFound exceptions and continues.

    :param monkeypatch: Pytest monkeypatch fixture
    :param patched_bq_client: Mock BigQuery client

    :raise: AssertionError

    :return: None
    """
    suffix = "staging_20250101000000_deadbe"
    specs = ingest_data_to_bigquery.define_ingestion_specs(suffix)

    def _maybe_raise(table_id: str) -> None:
        # Always record the attempt, raise NotFound on the first call only
        patched_bq_client.delete_table_recorder.append(table_id)
        if len(patched_bq_client.delete_table_recorder) == 1:
            raise gapi_exceptions.NotFound("missing")

    monkeypatch.setattr(patched_bq_client, "delete_table", _maybe_raise)

    ingest_data_to_bigquery.cleanup_staging_tables(
        client=patched_bq_client,
        project_id="proj",
        dataset="ds",
        ingestion_specs=specs,
    )

    # One failure + 3 successes → still recorded 4 attempts (the first raised but handled)
    assert len(patched_bq_client.delete_table_recorder) == 4
