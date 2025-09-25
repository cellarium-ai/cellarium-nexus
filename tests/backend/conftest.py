import pytest


@pytest.fixture(autouse=True)
def _enable_db_access(db):
    """
    Enable database access for cell_management unit tests.

    Database-backed cache (DatabaseCache) operations require an active
    database connection even in unit tests. Requesting the ``db`` fixture
    ensures Django test DB is set up.
    """
    pass
