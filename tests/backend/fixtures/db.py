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


@pytest.fixture(autouse=True)
def backend_test_db_setup(request: pytest.FixtureRequest):
    """
    Configure a file-backed SQLite database for backend tests.

    Use a real file so multiple connections share the same schema and data
    within the test session. The backend's base settings already configure
    DatabaseCache pointing at the "core_nexuscache" table; ensure that table
    exists for tests.
    """
    # Only run for tests under tests/backend/
    fspath = str(getattr(request.node, "fspath", ""))
    if "/tests/backend/" not in fspath:
        return

    from django.core.management import call_command
    from django.db import connections

    try:
        connections["default"].close()
    except Exception:
        pass
    call_command("createcachetable", "core_nexuscache", verbosity=0)
