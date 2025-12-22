import typing

import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import prepare_extract as prepare_extract_module
from cellarium.nexus.shared import schemas as schemas_module


class FakePreparer:
    """
    Provide a fake ExtractTablePreparer to record calls and arguments.
    """

    def __init__(self, *, client: typing.Any, project: str, dataset: str, extract_table_prefix: str) -> None:
        self.client = client
        self.project = project
        self.dataset = dataset
        self.prefix = extract_table_prefix
        self.calls: dict[str, list[tuple[tuple, dict]]] = {
            "prepare_feature_table": [],
            "prepare_cell_info": [],
            "prepare_count_matrix": [],
        }

    def prepare_feature_table(self, features: typing.Sequence[schemas_module.FeatureSchema]) -> None:  # noqa: D401
        self.calls["prepare_feature_table"].append(((features,), {}))

    def prepare_cell_info(
        self,
        *,
        extract_bin_size: int | None,
        random_seed_offset: int,
        partition_bin_count: int,
        partition_size: int,
        extract_bin_keys: list[str] | None,
        filters: dict[str, typing.Any] | None,
        obs_columns: list[str] | None,
        metadata_extra_columns: list[str] | None,
    ) -> None:
        self.calls["prepare_cell_info"].append(
            (
                tuple(),
                {
                    "extract_bin_size": extract_bin_size,
                    "random_seed_offset": random_seed_offset,
                    "partition_bin_count": partition_bin_count,
                    "partition_size": partition_size,
                    "extract_bin_keys": extract_bin_keys,
                    "filters": filters,
                    "obs_columns": obs_columns,
                    "metadata_extra_columns": metadata_extra_columns,
                },
            )
        )

    def prepare_count_matrix(self, *, partition_bin_count: int, partition_size: int) -> None:  # noqa: D401
        self.calls["prepare_count_matrix"].append(
            (tuple(), {"partition_bin_count": partition_bin_count, "partition_size": partition_size})
        )


class FakeMetadataExtractor:
    """
    Provide a fake MetadataExtractor to return deterministic compose_extract_metadata output.
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
        self.init_args = {
            "client": client,
            "project": project,
            "dataset": dataset,
            "extract_table_prefix": extract_table_prefix,
            "filters": filters,
            "extract_bin_size": extract_bin_size,
            "categorical_column_count_limit": categorical_column_count_limit,
        }

    def compose_extract_metadata(self) -> dict[str, typing.Any]:
        return {"ok": True, "source": "_FakeMetadataExtractor"}


def test_prepare_extract_tables_orchestrates_and_returns_metadata(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Validate that prepare_extract_tables orchestrates steps and returns metadata.

    Ensure that ExtractTablePreparer methods are invoked with correct arguments
    and that the final return value matches MetadataExtractor.compose_extract_metadata().

    :param monkeypatch: Pytest monkeypatch fixture
    """
    # Patch symbols inside module under test
    fake = FakePreparer(client=object(), project="proj", dataset="ds", extract_table_prefix="pref_")
    # Patch the class to return our prepared instance so we can inspect calls
    monkeypatch.setattr(prepare_extract_module, "ExtractTablePreparer", lambda **kwargs: fake)
    monkeypatch.setattr(prepare_extract_module, "MetadataExtractor", FakeMetadataExtractor)

    features = [
        schemas_module.FeatureSchema(id=1, symbol="X", ensemble_id="e1"),
        schemas_module.FeatureSchema(id=2, symbol="Y", ensemble_id="e2"),
    ]

    result = prepare_extract_module.prepare_extract_tables(
        client=object(),
        project="proj",
        dataset="ds",
        extract_table_prefix="pref_",
        features=features,
        categorical_column_count_limit=77,
        extract_bin_size=100,
        random_seed_offset=3,
        partition_bin_count=111,
        partition_size=8,
        extract_bin_keys=["k1"],
        filters={"a": 1},
        obs_columns=["c1", "c2"],
        metadata_extra_columns=["m1"],
    )

    # Validate feature table call
    assert len(fake.calls["prepare_feature_table"]) == 1
    (args, kwargs) = fake.calls["prepare_feature_table"][0]
    assert list(args[0]) == features  # same objects
    assert kwargs == {}

    # Validate cell info call wiring
    assert len(fake.calls["prepare_cell_info"]) == 1
    (_args2, kwargs2) = fake.calls["prepare_cell_info"][0]
    assert kwargs2 == {
        "extract_bin_size": 100,
        "random_seed_offset": 3,
        "partition_bin_count": 111,
        "partition_size": 8,
        "extract_bin_keys": ["k1"],
        "filters": {"a": 1},
        "obs_columns": ["c1", "c2"],
        "metadata_extra_columns": ["m1"],
    }

    # Validate count matrix call wiring
    assert len(fake.calls["prepare_count_matrix"]) == 1
    (_args3, kwargs3) = fake.calls["prepare_count_matrix"][0]
    assert kwargs3 == {"partition_bin_count": 111, "partition_size": 8}

    # Validate returned metadata
    assert result == {"ok": True, "source": "_FakeMetadataExtractor"}
