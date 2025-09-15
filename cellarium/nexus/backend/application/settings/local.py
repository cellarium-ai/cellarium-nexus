# flake8: noqa
from .base import *  # noqa

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
        },
    },
    "loggers": {
        "django.db.backends": {
            "handlers": ["console"],
            "level": "DEBUG",
        },
    },
}

# --- Local static files overrides ---
# In local development, serve static files directly from the app instead of GCS.
# This overrides STATIC_URL and the staticfiles storage backend to the built-in
# Django StaticFilesStorage, so changes to JS/CSS are picked up immediately.

# Use the local static URL
STATIC_URL = "/static/"

# Ensure the staticfiles storage uses the local filesystem instead of GCS
STORAGES["staticfiles"] = {
    "BACKEND": "django.contrib.staticfiles.storage.StaticFilesStorage",
}

STATIC_ROOT = BASE_DIR / "conf" / "static"
