import types

import pytest

from cellarium.nexus.omics_datastore.bq_sql.validation.query_validators import bigquery_validator as bq_val_module

cloud = pytest.importorskip("google.cloud")


class OkClient:
    def query(self, query: str, job_config: object):  # noqa: D401 - stub
        # return a minimal object that looks like a QueryJob
        return types.SimpleNamespace(state="DONE", query=query, job_config=job_config)


class ErrClient:
    def __init__(self, exc: Exception) -> None:
        self._exc = exc

    def query(self, query: str, job_config: object):  # noqa: D401 - stub
        raise self._exc


def test_bq_validator_success(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Validate that BigQuerySQLSyntaxValidator passes when the client returns a successful dry run.
    """
    monkeypatch.setattr(cloud.bigquery, "Client", lambda: OkClient())
    # should not raise
    bq_val_module.BigQuerySQLSyntaxValidator.validate_syntax(sql_query="select 1")


def test_bq_validator_raises_with_parsed_message(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Verify that a GoogleCloudError is parsed and surfaced with the expected message format.
    """

    # Craft an error whose str() matches the parsing logic used in the validator
    class FakeGCError(cloud.exceptions.GoogleCloudError):  # type: ignore[misc]
        def __str__(self) -> str:  # noqa: D401 - stub
            return "some-prefix jobs?prettyPrint=false: Syntax error: Unexpected keyword\n\n" "Trailing details..."

    monkeypatch.setattr(cloud.bigquery, "Client", lambda: ErrClient(FakeGCError("boom")))

    with pytest.raises(Exception) as ei:
        bq_val_module.BigQuerySQLSyntaxValidator.validate_syntax(sql_query="select broken")

    # the parsed message should be the first line after the marker and before blank line
    assert "Syntax error: Unexpected keyword" in str(ei.value)
