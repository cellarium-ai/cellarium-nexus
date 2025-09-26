from pathlib import Path

import environ

env = environ.Env()

BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent.parent.parent


env.read_env(BASE_DIR / "conf/.env")
# If running tests, load the dedicated test env file to override values
if env("ENVIRONMENT") == "test":
    env.read_env(BASE_DIR / "tests" / "fixtures" / "test_env", override=True)

if env("ENVIRONMENT") == "development":
    from .development import *  # noqa
elif env("ENVIRONMENT") == "local":
    try:
        from .local import *  # noqa
    except ImportError:
        from .development import *  # noqa
elif env("ENVIRONMENT") == "test":
    from .test import *  # noqa
else:
    from .production import *  # noqa
