import pytest


@pytest.fixture(scope="session", autouse=True)
def backend_test_db_setup():
    """
    Configure a file-backed SQLite database for backend tests.

    Use a real file so multiple connections share the same schema and data
    within the test session. The backend's base settings already configure
    DatabaseCache pointing at the "core_nexuscache" table; ensure that table
    exists for tests.
    """
    from django.core.management import call_command
    from django.db import connections

    # Ensure the cache table exists (createcachetable is idempotent)
    # We need a real DB connection available to run this management command.
    # For SQLite, the connection will create the file on demand.
    # Close any existing connection to ensure a clean state before creating the table.
    try:
        connections["default"].close()
    except Exception:
        pass
    call_command("createcachetable", "core_nexuscache", verbosity=0)
