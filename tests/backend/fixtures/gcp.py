from __future__ import annotations

import pytest
from google.cloud import bigquery, storage


class DummyBigQueryClient:
    """
    Represent a no-op BigQuery client used in backend tests.
    """

    def __init__(self, *args, **kwargs) -> None:  # noqa: D401
        """Accept any arguments without side effects."""


@pytest.fixture(autouse=True)
def dummy_bigquery_client(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Patch bigquery.Client to a no-op dummy for backend tests.

    This prevents network or ADC interactions when code constructs a client
    internally (for example, ``OmicsCachedDataManager`` with ``BigQueryDataOperator``).
    """

    monkeypatch.setattr(bigquery, "Client", DummyBigQueryClient, raising=True)


class DummyGCSClient:
    """
    Provide a minimal stub for Google Cloud Storage interactions.

    :param bucket_name: Optional bucket name to bind for path assertions
    """

    def __init__(self, bucket_name: str | None = None) -> None:
        self.bucket_name = bucket_name
        self.calls: list[tuple[str, str]] = []

    class _Bucket:
        def __init__(self, name: str, parent: "DummyGCSClient") -> None:
            self.name = name
            self._parent = parent

        def blob(self, path: str) -> "DummyGCSClient._Blob":
            return DummyGCSClient._Blob(bucket=self, path=path)

        def list_blobs(self, prefix: str | None = None, max_results: int | None = None):
            """
            Provide a minimal iterable of blobs. By default returns empty.

            :param prefix: Optional prefix filter (ignored in dummy)
            :param max_results: Optional maximum number of results (ignored in dummy)

            :return: Empty list to indicate no objects found
            """
            return []

    class _Blob:
        def __init__(self, bucket: "DummyGCSClient._Bucket", path: str) -> None:
            self.bucket = bucket
            self.path = path

        def upload_from_filename(self, filename: str, *_, **__) -> None:
            self.bucket._parent.calls.append((self.path, filename))

    def bucket(self, name: str) -> "DummyGCSClient._Bucket":
        return DummyGCSClient._Bucket(name=name, parent=self)


@pytest.fixture()
def dummy_gcs_client(monkeypatch: pytest.MonkeyPatch) -> DummyGCSClient:
    """
    Patch ``google.cloud.storage.Client`` with an in-memory stub.

    :param monkeypatch: Pytest monkeypatch fixture

    :return: Dummy GCS client instance recording uploads
    """

    client = DummyGCSClient()
    monkeypatch.setattr(storage, "Client", lambda *_, **__: client, raising=True)
    return client


class VertexAIPipelineStub:
    """
    Record submissions to Vertex AI pipeline helper.
    """

    def __init__(self) -> None:
        self.calls: list[dict[str, object]] = []
        self.default_url = "https://vertex.example/runs/stub"

    def __call__(self, *args: object, **kwargs: object) -> str:
        self.calls.append({"args": args, "kwargs": kwargs})
        return self.default_url


@pytest.fixture()
def vertex_ai_pipeline_stub(monkeypatch: pytest.MonkeyPatch) -> VertexAIPipelineStub:
    """
    Patch ``submit_pipeline`` helper to capture Vertex AI invocations.

    :param monkeypatch: Pytest monkeypatch fixture

    :return: Callable stub retaining submission metadata
    """

    stub = VertexAIPipelineStub()
    monkeypatch.setattr(
        "cellarium.nexus.workflows.kubeflow.utils.job.submit_pipeline",
        stub,
        raising=True,
    )
    return stub
