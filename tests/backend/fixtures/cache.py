import typing

import pytest


@pytest.fixture(autouse=True)
def backend_cache_cleaner(settings: typing.Any) -> None:
    """
    Clear the configured cache (DatabaseCache via SQLite) around each test for isolation.

    :param settings: Django settings fixture
    """
    # Clear cache between tests for isolation
    from django.core.cache import cache

    cache.clear()
    yield
    cache.clear()
