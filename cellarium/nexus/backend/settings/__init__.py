from pathlib import Path

import environ

env = environ.Env()

BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent.parent


env.read_env(BASE_DIR / "cellarium/nexus/.env")

if env("ENVIRONMENT") == "development":
    from .development import *  # noqa
elif env("ENVIRONMENT") == "local":
    try:
        from .local import *  # noqa
    except ImportError:
        from .development import *  # noqa
else:
    from .production import *  # noqa
