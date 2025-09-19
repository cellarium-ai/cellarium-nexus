import datetime as dt
import typing

import google.auth as google_auth
import google.auth.credentials as google_auth_credentials
import pytest
from google.cloud import bigquery


@pytest.fixture()
def freeze_time(monkeypatch: pytest.MonkeyPatch) -> typing.Callable[[dt.datetime], None]:
    """
    Freeze datetime.utcnow() and datetime.now(datetime.UTC) to a fixed value.

    Provide a setter to change frozen time during a test.

    :param monkeypatch: Pytest monkeypatch fixture

    :raise: None

    :return: Setter function that accepts a datetime value
    """
    frozen: dict[str, dt.datetime] = {"value": dt.datetime(2025, 1, 1, 0, 0, 0)}

    class _FrozenDateTime(dt.datetime):
        @classmethod
        def utcnow(cls) -> dt.datetime:  # type: ignore[override]
            return frozen["value"]

        @classmethod
        def now(cls, tz: dt.tzinfo | None = None) -> dt.datetime:  # type: ignore[override]
            val = frozen["value"]
            if tz is None:
                return val
            return val.replace(tzinfo=tz)

    monkeypatch.setattr(dt, "datetime", _FrozenDateTime)

    # Python 3.11+ introduces datetime.UTC; provide a compatibility shim for 3.10
    if not hasattr(dt, "UTC"):
        monkeypatch.setattr(dt, "UTC", dt.timezone.utc, raising=False)

    def _set(new_value: dt.datetime) -> None:
        frozen["value"] = new_value

    return _set


@pytest.fixture()
def freeze_uuid(monkeypatch: pytest.MonkeyPatch) -> typing.Callable[[str], None]:
    """
    Freeze uuid.uuid4().hex to a chosen 32-char hex string.

    :param monkeypatch: Pytest monkeypatch fixture

    :raise: None

    :return: Setter function that accepts a 32-hex string (lowercase)
    """
    import uuid as _uuid

    state: dict[str, str] = {"hex": "deadbeefcafebabe0123456789abcd0"}

    class _FrozenUUID:
        def __init__(self, hex: str) -> None:
            self.hex = hex

    def _uuid4() -> _FrozenUUID:
        return _FrozenUUID(state["hex"])

    monkeypatch.setattr(_uuid, "uuid4", _uuid4)

    def _set(hex32: str) -> None:
        state["hex"] = hex32

    return _set


class BQLoadJobMock:
    """
    Represent a dummy BigQuery Load Job with minimal surface for tests.

    Provide ``result()`` and ``output_rows`` attributes similar to real job.
    """

    def __init__(self) -> None:
        self.output_rows: int = 0

    def result(self) -> None:
        return None


class BQQueryJobMock:
    """
    Represent a dummy BigQuery Query Job with minimal surface for tests.

    Provide ``result()`` similar to real job.
    """

    def result(self) -> None:
        return None


class BQClientMock:
    """
    Represent a mock BigQuery client with inspectable behavior for tests.

    Attributes adjustable in tests:
    - load_job_rows: int
    - query_sql_recorder: list[str]
    - delete_table_recorder: list[str]
    - last_load_uri: str | None
    - last_load_table_id: str | None
    - last_load_job_config: bigquery.LoadJobConfig | None
    """

    def __init__(self) -> None:
        self.load_job_rows: int = 0
        self.query_sql_recorder: list[str] = []
        self.delete_table_recorder: list[str] = []
        self.last_load_uri: str | None = None
        self.last_load_table_id: str | None = None
        self.last_load_job_config: bigquery.LoadJobConfig | None = None
        # Job instances injected by fixture
        self._bq_load_job: BQLoadJobMock | None = None
        self._bq_query_job: BQQueryJobMock | None = None

    def load_table_from_uri(
        self,
        uri: str,
        table_id: str,
        job_config: bigquery.LoadJobConfig,
    ) -> BQLoadJobMock:
        """
        Record parameters and return a provided load job mock.

        :param uri: GCS URI of the source data
        :param table_id: Target BigQuery table ID
        :param job_config: BigQuery LoadJobConfig

        :raise: None

        :return: Provided load job mock instance
        """
        self.last_load_uri = uri
        self.last_load_table_id = table_id
        self.last_load_job_config = job_config
        job = self._bq_load_job or BQLoadJobMock()
        # update job instance rows for assertion
        setattr(job, "output_rows", self.load_job_rows)
        return job

    def query(self, sql: str) -> BQQueryJobMock:
        """
        Record SQL and return the provided query job mock.

        :param sql: SQL text to record

        :raise: None

        :return: Provided query job mock instance
        """
        self.query_sql_recorder.append(sql)
        return self._bq_query_job or BQQueryJobMock()

    def delete_table(self, table_id: str) -> None:
        """
        Record deleted table ID for assertions.

        :param table_id: Table to delete

        :raise: None
        """
        self.delete_table_recorder.append(table_id)


@pytest.fixture(autouse=True)
def disable_google_adc(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Disable Google Application Default Credentials lookup in tests.

    Patch ``google.auth.default`` to return anonymous credentials and a dummy
    project to prevent network calls and credential discovery in CI.

    :param monkeypatch: Pytest monkeypatch fixture

    :raise: None

    :return: None
    """
    creds = google_auth_credentials.AnonymousCredentials()

    def _default(*args: typing.Any, **kwargs: typing.Any) -> tuple[google_auth_credentials.AnonymousCredentials, str]:
        return creds, "test-project"

    monkeypatch.setattr(target=google_auth, name="default", value=_default)


@pytest.fixture()
def bq_load_job() -> BQLoadJobMock:
    """
    Provide a dummy BigQuery Load Job instance with result() and output_rows.

    :raise: None

    :return: Load job instance
    """
    return BQLoadJobMock()


@pytest.fixture()
def bq_query_job() -> BQQueryJobMock:
    """
    Provide a dummy BigQuery Query Job instance with result().

    :raise: None

    :return: Query job instance
    """
    return BQQueryJobMock()


@pytest.fixture()
def bq_client(bq_load_job: BQLoadJobMock, bq_query_job: BQQueryJobMock) -> BQClientMock:
    """
    Provide a mock BigQuery client instance with simple, inspectable behavior.

    Instance attributes you can adjust in tests:
    - load_job_rows: int
    - query_sql_recorder: list[str]
    - delete_table_recorder: list[str]
    - last_load_uri: str | None
    - last_load_table_id: str | None
    - last_load_job_config: bigquery.LoadJobConfig | None

    :raise: None

    :return: Mock client instance
    """
    client = BQClientMock()
    # Inject job mocks for use by methods
    client._bq_load_job = bq_load_job
    client._bq_query_job = bq_query_job
    return client
