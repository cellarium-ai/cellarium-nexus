# Configure Django settings for tests in case pytest-django plugin is not installed/loaded.
# This must run before importing/activating any test plugins that may import Django code.

pytest_plugins = [
    "tests.fixtures.data_inputs",
    "tests.fixtures.helpers",
    "tests.fixtures.mocks",
]
