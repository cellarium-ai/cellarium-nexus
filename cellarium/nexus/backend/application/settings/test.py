# ruff: noqa
from .base import *  # noqa

# Testing overrides
DEBUG = False

# Use a file-backed SQLite DB for tests to support multi-connection behavior
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": str(BASE_DIR / "conf" / "test.sqlite3"),
    }
}


LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "stream": "ext://sys.stdout",
        },
    },
    "loggers": {
        "django.db.backends": {
            "handlers": ["console"],
            "level": "DEBUG",  # SQL statements
            "propagate": False,
        },
        # Django internals
        "django": {"handlers": ["console"], "level": "INFO", "propagate": False},
    },
}
