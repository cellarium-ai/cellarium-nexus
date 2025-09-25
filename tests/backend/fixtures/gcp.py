import pytest
from google.cloud import bigquery


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
    internally (for example, ``BigQueryCachedDataManager``).
    """

    monkeypatch.setattr(bigquery, "Client", DummyBigQueryClient, raising=True)
