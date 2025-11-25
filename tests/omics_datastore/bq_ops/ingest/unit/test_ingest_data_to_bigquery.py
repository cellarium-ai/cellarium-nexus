import pathlib
import re
from typing import Any

import pytest
from google.api_core import exceptions as gapi_exceptions
from google.cloud import bigquery

from cellarium.nexus.omics_datastore.bq_ops import constants, exceptions
from cellarium.nexus.omics_datastore.bq_ops.ingest import ingest_data_to_bigquery as bq_ingest


def test_initialize_bigquery_client(monkeypatch: pytest.MonkeyPatch, bq_client: Any) -> None:
    """
    Initialize and return a BigQuery client using provided project id.
    """
    created: dict[str, Any] = {}

    def _make_client(project: str) -> Any:
        created["project"] = project
        return bq_client

    monkeypatch.setattr(bq_ingest.bigquery, "Client", _make_client)

    client = bq_ingest.initialize_bigquery_client(project_id="proj1", dataset="ds")
    assert client is bq_client
    assert created == {"project": "proj1"}


def test_generate_staging_suffix() -> None:
    """
    Generate a unique staging suffix matching the expected pattern.
    """
    suffix = bq_ingest.generate_staging_suffix()
    assert re.match(r"^staging_\d{14}_[0-9a-f]{6}$", suffix)


def test_define_ingestion_specs() -> None:
    """
    Define ingestion specs for all expected tables with correct formats and names.
    """
    suffix = "staging_20240101123456_ab12cd"
    specs = bq_ingest.define_ingestion_specs(staging_suffix=suffix)

    assert len(specs) == 4
    # Check table names and formats
    expected = [
        (constants.BQ_INGEST_TABLE_NAME, bigquery.SourceFormat.AVRO),
        (constants.BQ_CELL_INFO_TABLE_NAME, bigquery.SourceFormat.AVRO),
        (constants.BQ_FEATURE_INFO_TABLE_NAME, bigquery.SourceFormat.AVRO),
        (constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME, bigquery.SourceFormat.CSV),
    ]
    for (final_name, _staging, pattern, fmt), (exp_name, exp_fmt) in zip(specs, expected):
        assert final_name == exp_name
        assert fmt == exp_fmt
        # File pattern matches constants used by implementation
        assert pattern in {
            constants.INGEST_INGEST_FILE_NAME,
            constants.INGEST_CELL_INFO_FILE_NAME,
            constants.INGEST_FEATURE_INFO_FILE_NAME,
            constants.INGEST_RAW_COUNTS_FILE_PATTERN,
        }


def test_load_table_from_gcs_happy_path(bq_client: Any) -> None:
    """
    Load table from GCS constructs correct URI/table_id and waits for job completion.
    """
    # Configure expected rows and call
    bq_client.load_job_rows = 42

    bq_ingest.load_table_from_gcs(
        client=bq_client,  # type: ignore[arg-type]
        project="proj",
        dataset="ds",
        table_name="t",
        file_pattern="file.avro",
        file_format=bigquery.SourceFormat.AVRO,
        gcs_bucket_name="bucket",
        gcs_stage_dir="stage",
    )

    assert bq_client.last_load_uri == "gs://bucket/stage/file.avro"
    assert bq_client.last_load_table_id == "proj.ds.t"
    assert isinstance(bq_client.last_load_job_config, bigquery.LoadJobConfig)


def test_perform_load_table_from_gcs_delegates(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Perform load calls underlying load function with same parameters.
    """
    calls: list[dict[str, Any]] = []

    def _fake_loader(**kwargs: Any) -> None:
        calls.append(kwargs)

    monkeypatch.setattr(bq_ingest, "load_table_from_gcs", _fake_loader)

    bq_ingest.perform_load_table_from_gcs(
        client=None,  # type: ignore[arg-type]
        project="p",
        dataset="d",
        table_name="t",
        file_pattern="*.avro",
        file_format=bigquery.SourceFormat.AVRO,
        gcs_bucket_name="b",
        gcs_stage_dir="s",
    )

    assert len(calls) == 1
    assert calls[0]["project"] == "p"
    assert calls[0]["dataset"] == "d"
    assert calls[0]["table_name"] == "t"


def test_load_data_into_staging_all_success(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Return True when all staging loads succeed.
    """
    called: list[tuple[str, str]] = []

    def _ok(**kwargs: Any) -> None:
        called.append((kwargs["table_name"], kwargs["file_pattern"]))

    monkeypatch.setattr(bq_ingest, "perform_load_table_from_gcs", _ok)

    specs = bq_ingest.define_ingestion_specs(staging_suffix="staging_test_abc123")
    ok = bq_ingest.load_data_into_staging(
        client=None,  # type: ignore[arg-type]
        project_id="p",
        dataset="d",
        gcs_bucket_name="b",
        gcs_stage_dir="s",
        ingestion_specs=specs,
    )
    assert ok is True
    assert len(called) == len(specs)


def test_load_data_into_staging_partial_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Return False when any staging load raises an exception.
    """

    def _sometimes_fail(**kwargs: Any) -> None:
        if kwargs["table_name"].endswith("feature_info_staging"):
            raise RuntimeError("boom")

    monkeypatch.setattr(bq_ingest, "perform_load_table_from_gcs", _sometimes_fail)

    specs = [
        (constants.BQ_INGEST_TABLE_NAME, "ingest_staging", "ingest-info.avro", bigquery.SourceFormat.AVRO),
        (constants.BQ_FEATURE_INFO_TABLE_NAME, "feature_info_staging", "feature-info.avro", bigquery.SourceFormat.AVRO),
    ]

    ok = bq_ingest.load_data_into_staging(
        client=None,  # type: ignore[arg-type]
        project_id="p",
        dataset="d",
        gcs_bucket_name="b",
        gcs_stage_dir="s",
        ingestion_specs=specs,
    )
    assert ok is False


def test_commit_to_production_success(bq_client: Any) -> None:
    """
    Build and execute multi-statement transaction and return True on success.
    """
    specs = [
        (
            constants.BQ_INGEST_TABLE_NAME,
            "ingest_staging",
            constants.INGEST_INGEST_FILE_NAME,
            bigquery.SourceFormat.AVRO,
        ),
        (
            constants.BQ_CELL_INFO_TABLE_NAME,
            "cell_info_staging",
            constants.INGEST_CELL_INFO_FILE_NAME,
            bigquery.SourceFormat.AVRO,
        ),
    ]
    ok = bq_ingest.commit_to_production(client=bq_client, dataset="ds", ingestion_specs=specs)
    assert ok is True
    assert any("BEGIN TRANSACTION" in sql for sql in bq_client.query_sql_recorder)
    assert any("COMMIT TRANSACTION" in sql for sql in bq_client.query_sql_recorder)
    assert any("INSERT INTO`ds.ingest".replace(" ", "") in sql.replace(" ", "") for sql in bq_client.query_sql_recorder)


def test_commit_to_production_failure(monkeypatch: pytest.MonkeyPatch, bq_client: Any) -> None:
    """
    Return False if query fails.
    """

    def _raise_query(sql: str) -> Any:
        raise RuntimeError("query error")

    monkeypatch.setattr(bq_client, "query", _raise_query)

    ok = bq_ingest.commit_to_production(client=bq_client, dataset="ds", ingestion_specs=[])
    assert ok is False


def test_cleanup_staging_tables(monkeypatch: pytest.MonkeyPatch, bq_client: Any) -> None:
    """
    Attempt to delete each staging table; ignore NotFound and log other errors.
    """
    specs = [
        ("final1", "staging1", "pat", bigquery.SourceFormat.AVRO),
        ("final2", "staging_miss", "pat", bigquery.SourceFormat.AVRO),
        ("final3", "staging_err", "pat", bigquery.SourceFormat.AVRO),
    ]

    def _delete_table(table_id: str) -> None:
        if table_id.endswith("miss"):
            raise gapi_exceptions.NotFound("missing")
        if table_id.endswith("err"):
            raise RuntimeError("delete error")
        bq_client.delete_table_recorder.append(table_id)

    monkeypatch.setattr(bq_client, "delete_table", _delete_table)

    bq_ingest.cleanup_staging_tables(client=bq_client, project_id="p", dataset="d", ingestion_specs=specs)
    assert bq_client.delete_table_recorder == ["p.d.staging1"]


def test_get_schema_for_table_success(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Fetch BigQuery schema from pydantic model via converter for a known table name.
    """
    sentinel_schema: list[bigquery.SchemaField] = [bigquery.SchemaField("c", "STRING")]

    def _fake_convert(pydantic_model: Any) -> list[bigquery.SchemaField]:
        return sentinel_schema

    monkeypatch.setattr(bq_ingest.converter, "pydantic_to_bigquery", _fake_convert)

    out = bq_ingest.get_schema_for_table(final_table=constants.BQ_CELL_INFO_TABLE_NAME)
    assert out is sentinel_schema


def test_get_schema_for_table_error() -> None:
    """
    Raise when table name is unknown.
    """
    with pytest.raises(exceptions.BigQuerySchemaError):
        _ = bq_ingest.get_schema_for_table(final_table="unknown")


def test_bigquery_ingest_context_happy_path(
    monkeypatch: pytest.MonkeyPatch, bq_client: Any, tmp_path: pathlib.Path
) -> None:
    """
    Orchestrate staging create+load then commit when body yields without error and cleanup finally.
    """
    # Fixed suffix for determinism
    monkeypatch.setattr(bq_ingest, "generate_staging_suffix", lambda: "staging_fixed_suffix")

    monkeypatch.setattr(bq_ingest, "initialize_bigquery_client", lambda project_id, dataset: bq_client)

    # Stub create staging
    created: list[tuple[str, str]] = []

    def _create_staging_table(
        client: Any,
        project: str,
        dataset: str,
        base_table_name: str,
        staging_suffix: str,
        schema: list[Any],
        clustering_fields: list[str] | None,
    ) -> None:  # noqa: E501
        created.append((base_table_name, staging_suffix))

    monkeypatch.setattr(bq_ingest, "create_staging_table", _create_staging_table)

    # Stub schema conversion
    monkeypatch.setattr(bq_ingest, "get_schema_for_table", lambda final_table: [])

    # Loads succeed
    monkeypatch.setattr(bq_ingest, "load_data_into_staging", lambda **_: True)

    # Commit succeeds
    monkeypatch.setattr(bq_ingest, "commit_to_production", lambda **_: True)

    with bq_ingest.bigquery_ingest_context(
        project_id="p",
        dataset="d",
        gcs_bucket_name="b",
        gcs_stage_dir="s",
    ):
        # Do nothing inside
        pass

    # Verify staging tables created for all final tables in specs
    assert any(t[0] == constants.BQ_INGEST_TABLE_NAME for t in created)
    assert any(t[0] == constants.BQ_CELL_INFO_TABLE_NAME for t in created)


def test_bigquery_ingest_context_load_failure(monkeypatch: pytest.MonkeyPatch, bq_client: Any) -> None:
    """
    Raise BigQueryLoadError when staging load reports failure; still cleanup.
    """
    monkeypatch.setattr(bq_ingest, "initialize_bigquery_client", lambda project_id, dataset: bq_client)
    monkeypatch.setattr(bq_ingest, "generate_staging_suffix", lambda: "sfx")
    monkeypatch.setattr(bq_ingest, "get_schema_for_table", lambda final_table: [])
    monkeypatch.setattr(bq_ingest, "create_staging_table", lambda **_: None)
    monkeypatch.setattr(bq_ingest, "load_data_into_staging", lambda **_: False)

    with pytest.raises(exceptions.BigQueryLoadError):
        with bq_ingest.bigquery_ingest_context(
            project_id="p",
            dataset="d",
            gcs_bucket_name="b",
            gcs_stage_dir="s",
        ):
            pass


def test_bigquery_ingest_context_commit_failure(monkeypatch: pytest.MonkeyPatch, bq_client: Any) -> None:
    """
    Raise BigQueryCommitError when commit returns False; still cleanup.
    """
    monkeypatch.setattr(bq_ingest, "initialize_bigquery_client", lambda project_id, dataset: bq_client)
    monkeypatch.setattr(bq_ingest, "generate_staging_suffix", lambda: "sfx")
    monkeypatch.setattr(bq_ingest, "get_schema_for_table", lambda final_table: [])
    monkeypatch.setattr(bq_ingest, "create_staging_table", lambda **_: None)
    monkeypatch.setattr(bq_ingest, "load_data_into_staging", lambda **_: True)
    monkeypatch.setattr(bq_ingest, "commit_to_production", lambda **_: False)

    with pytest.raises(exceptions.BigQueryCommitError):
        with bq_ingest.bigquery_ingest_context(
            project_id="p",
            dataset="d",
            gcs_bucket_name="b",
            gcs_stage_dir="s",
        ):
            pass
