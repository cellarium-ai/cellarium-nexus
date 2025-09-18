import typing

import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import prepare_extract as prepare_extract_module
from cellarium.nexus.shared import schemas as schemas_module


class FakeLoadJob:
    """
    Represent a dummy load job with result() and output_rows attributes.
    """

    def __init__(self, output_rows: int = 2) -> None:
        self.output_rows = output_rows

    def result(self) -> None:
        return None


class Tbl:
    """
    Minimal table stub with an 'expires' attribute used by update_table.
    """

    def __init__(self) -> None:
        self.expires = None  # type: ignore[assignment]


class FakeMetadataExtractor:
    """
    Deterministic metadata extractor stub for integration-style test.
    """

    def __init__(
        self,
        *,
        client: typing.Any,
        project: str,
        dataset: str,
        extract_table_prefix: str,
        filters: dict[str, typing.Any] | None,
        extract_bin_size: int | None,
        categorical_column_count_limit: int,
    ) -> None:
        self._args = {
            "client": client,
            "project": project,
            "dataset": dataset,
            "extract_table_prefix": extract_table_prefix,
            "filters": filters,
            "extract_bin_size": extract_bin_size,
            "categorical_column_count_limit": categorical_column_count_limit,
        }

    def compose_extract_metadata(self) -> dict[str, typing.Any]:
        return {"ok": True}


def test_prepare_extract_tables_happy_path(
    monkeypatch: pytest.MonkeyPatch,
    bq_client: typing.Any,
) -> None:
    """
    Exercise prepare_extract_tables end-to-end with a BigQuery client mock.

    Validate that:
    - Feature table data is loaded via load_table_from_file
    - Three SQL queries run for randomized, binned, and drop steps
    - Table expiration updates occur for feature, cell info, and count matrix tables
    - Returned metadata comes from MetadataExtractor.compose_extract_metadata

    :param monkeypatch: Pytest monkeypatch fixture
    :param bq_client: BigQuery client mock fixture
    """
    # Capture bq_sql.render calls and return recognizable SQL
    render_calls: list[tuple[str, typing.Any]] = []

    def _render(path: str, data: typing.Any) -> str:
        render_calls.append((path, data))
        return f"-- {path}\nSELECT 1;"

    monkeypatch.setattr(prepare_extract_module.bq_sql, "render", _render)

    # Add methods missing on BQClientMock
    bq_client._ltf_calls: list[tuple[str, typing.Any]] = []
    bq_client._update_calls: list[list[str]] = []

    def _load_table_from_file(file_obj: typing.Any, table_id: str, job_config: typing.Any) -> typing.Any:
        bq_client._ltf_calls.append((table_id, job_config))
        return FakeLoadJob()

    def _get_table(table_id: str) -> Tbl:
        return Tbl()

    def _update_table(table: Tbl, fields: list[str]) -> None:
        bq_client._update_calls.append(fields)

    monkeypatch.setattr(bq_client, "load_table_from_file", _load_table_from_file, raising=False)
    monkeypatch.setattr(bq_client, "get_table", _get_table, raising=False)
    monkeypatch.setattr(bq_client, "update_table", _update_table, raising=False)

    # Stub MetadataExtractor to deterministic output
    monkeypatch.setattr(prepare_extract_module, "MetadataExtractor", FakeMetadataExtractor)

    features = [
        schemas_module.FeatureSchema(id=1, symbol="A", ensemble_id="e1"),
        schemas_module.FeatureSchema(id=2, symbol="B", ensemble_id="e2"),
    ]

    result = prepare_extract_module.prepare_extract_tables(
        client=bq_client,
        project="proj",
        dataset="ds",
        extract_table_prefix="pref_",
        features=features,
        categorical_column_count_limit=50,
        extract_bin_size=100,
        random_seed_offset=5,
        partition_bin_count=10,
        partition_size=2,
        extract_bin_keys=["k1", "k2"],
        filters={"organism__eq": "homo"},
        obs_columns=["c1", "c2"],
        metadata_extra_columns=["m1"],
    )

    # 1) Feature load issued
    assert len(bq_client._ltf_calls) == 1
    table_id, job_config = bq_client._ltf_calls[0]
    assert table_id.startswith("proj.ds.pref_")
    assert getattr(job_config, "source_format", None) is not None

    # 2) Four queries executed (rand, binned, drop, count_matrix)
    assert len(bq_client.query_sql_recorder) == 4
    # Also ensure render was invoked for all 4 templates (two cell_info, one drop, one count_matrix)
    assert len(render_calls) == 4

    # 3) Expiration updates for feature, cell info, count matrix (3 updates)
    assert len(bq_client._update_calls) == 3
    assert all(fields == ["expires"] for fields in bq_client._update_calls)

    # 4) Metadata returned
    assert result == {"ok": True}
