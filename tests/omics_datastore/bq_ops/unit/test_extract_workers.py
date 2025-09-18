import concurrent.futures as cf
import typing
from pathlib import Path

import pytest

from cellarium.nexus.omics_datastore.bq_ops.extract import extract as extract_module
from cellarium.nexus.shared import schemas as schemas_module


class FakeExtractorForPerform:
    """
    Record calls to extract_bin_to_anndata for perform_extraction tests.
    """

    def __init__(self) -> None:
        self.calls: list[dict[str, typing.Any]] = []

    def extract_bin_to_anndata(
        self,
        *,
        bin_number: int,
        output_path: Path,
        extract_metadata: schemas_module.ExtractMetadata,
        obs_columns: list[str] | None = None,
    ) -> None:
        self.calls.append(
            {
                "bin_number": bin_number,
                "output_path": output_path,
                "extract_metadata": extract_metadata,
                "obs_columns": obs_columns,
            }
        )


class FakeClientTracking:
    """
    Track BigQuery client construction and close calls.
    """

    constructed_project: str | None = None
    closed_count: int = 0

    def __init__(self, project: str) -> None:
        FakeClientTracking.constructed_project = project

    def close(self) -> None:  # noqa: D401
        FakeClientTracking.closed_count += 1


class FakeExtractorForWorker:
    """
    Track DataExtractor construction args and performed extraction kwargs in worker tests.
    """

    constructed: dict[str, typing.Any] = {}
    performed: dict[str, typing.Any] = {}

    def __init__(self, *, client: typing.Any, project: str, dataset: str, extract_table_prefix: str) -> None:
        FakeExtractorForWorker.constructed = {"dataset": dataset, "prefix": extract_table_prefix}

    def extract_bin_to_anndata(self, **kwargs: typing.Any) -> None:  # noqa: D401
        FakeExtractorForWorker.performed = kwargs


class FakeExecutorOK:
    """
    Provide an executor whose futures complete successfully.
    """

    def __init__(self, *args: typing.Any, **kwargs: typing.Any) -> None:
        pass

    def __enter__(self) -> "FakeExecutorOK":  # noqa: D401
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # noqa: D401
        return None

    def submit(self, fn: typing.Callable[..., None], **kwargs: typing.Any):  # noqa: D401
        FakeExecutorOK.submitted.append(kwargs)
        fut: cf.Future[None] = cf.Future()
        fut.set_result(None)
        return fut

    # store submissions for assertions
    submitted: list[dict[str, typing.Any]] = []


class FakeExecutorError:
    """
    Provide an executor whose futures raise an error when awaited.
    """

    def __init__(self, *args: typing.Any, **kwargs: typing.Any) -> None:
        pass

    def __enter__(self) -> "FakeExecutorError":  # noqa: D401
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # noqa: D401
        return None

    def submit(self, fn: typing.Callable[..., None], **kwargs: typing.Any):  # noqa: D401
        fut: cf.Future[None] = cf.Future()
        fut.set_exception(RuntimeError("boom"))
        return fut


def test_perform_extraction_calls_underlying_once(tmp_path: Path) -> None:
    """
    Verify that perform_extraction invokes extractor.extract_bin_to_anndata once.

    :param tmp_path: Temporary directory path provided by pytest
    """
    extractor = FakeExtractorForPerform()

    meta = schemas_module.ExtractMetadata(total_bins=1, last_bin_size=1)
    extract_module.perform_extraction(
        extractor=extractor,
        bin_number=7,
        output_path=tmp_path / "x.h5ad",
        extract_metadata=meta,
        obs_columns=["c1"],
    )

    assert len(extractor.calls) == 1
    assert extractor.calls[0]["bin_number"] == 7


def test_extract_bin_to_anndata_worker_constructs_client_and_closes(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """
    Ensure extract_bin_to_anndata_worker constructs BigQuery client, performs extraction, and closes client.

    :param monkeypatch: Pytest monkeypatch fixture
    :param tmp_path: Temporary directory path provided by pytest
    """
    # reset trackers
    FakeClientTracking.constructed_project = None
    FakeClientTracking.closed_count = 0
    FakeExtractorForWorker.constructed = {}
    FakeExtractorForWorker.performed = {}

    monkeypatch.setattr(extract_module.bigquery, "Client", FakeClientTracking)
    monkeypatch.setattr(extract_module, "DataExtractor", FakeExtractorForWorker)

    meta = schemas_module.ExtractMetadata(total_bins=1, last_bin_size=1)

    extract_module.extract_bin_to_anndata_worker(
        project="P",
        dataset="D",
        extract_table_prefix="X_",
        bin_number=9,
        output_path=tmp_path / "out.h5ad",
        extract_metadata=meta,
        obs_columns=None,
    )

    assert FakeClientTracking.constructed_project == "P"
    assert FakeExtractorForWorker.constructed["dataset"] == "D"
    assert FakeExtractorForWorker.constructed["prefix"] == "X_"
    assert FakeClientTracking.closed_count == 1
    assert FakeExtractorForWorker.performed["bin_number"] == 9


def test_extract_bins_submits_jobs_and_uses_output_names(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    patched_process_pool_executor: typing.Iterator[None],
) -> None:
    """
    Check that extract_bins submits one job per bin and writes expected output file names.

    :param monkeypatch: Pytest monkeypatch fixture
    :param tmp_path: Temporary directory path provided by pytest
    :param patched_process_pool_executor: Fixture patching ProcessPoolExecutor to threads
    """
    # reset submissions
    FakeExecutorOK.submitted = []
    monkeypatch.setattr(extract_module.concurrency, "ProcessPoolExecutor", FakeExecutorOK)

    meta = schemas_module.ExtractMetadata(total_bins=2, last_bin_size=1)

    extract_module.extract_bins(
        client=object(),
        project="P",
        dataset="D",
        extract_table_prefix="X_",
        bins=[1, 2],
        output_dir=tmp_path,
        extract_metadata=meta,
        obs_columns=None,
        max_workers=2,
    )

    submitted = FakeExecutorOK.submitted
    assert len(submitted) == 2
    assert submitted[0]["bin_number"] == 1
    assert submitted[0]["output_path"].name == "extract_1.h5ad"
    assert submitted[1]["bin_number"] == 2
    assert submitted[1]["output_path"].name == "extract_2.h5ad"


def test_extract_bins_raises_on_worker_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    patched_process_pool_executor: typing.Iterator[None],
) -> None:
    """
    Validate that extract_bins surfaces exceptions raised by worker futures.

    :param monkeypatch: Pytest monkeypatch fixture
    :param tmp_path: Temporary directory path provided by pytest
    :param patched_process_pool_executor: Fixture patching ProcessPoolExecutor to threads
    """
    monkeypatch.setattr(extract_module.concurrency, "ProcessPoolExecutor", FakeExecutorError)

    meta = schemas_module.ExtractMetadata(total_bins=2, last_bin_size=1)

    with pytest.raises(RuntimeError):
        extract_module.extract_bins(
            client=object(),
            project="P",
            dataset="D",
            extract_table_prefix="X_",
            bins=[1, 2],
            output_dir=tmp_path,
            extract_metadata=meta,
            obs_columns=None,
            max_workers=2,
        )
