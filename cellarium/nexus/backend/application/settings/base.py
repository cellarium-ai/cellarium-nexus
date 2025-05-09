"""Base Django settings for Cellarium Nexus project."""

from django.templatetags.static import static
from django.urls import reverse_lazy
from django.utils.translation import gettext_lazy as _

from cellarium.nexus.backend.application.settings import BASE_DIR, env

# Core Django Settings
# -------------------
SECRET_KEY = env("SECRET_KEY")
DEBUG = True
ALLOWED_HOSTS = [env("MAIN_HOST_ALLOWED")]
SITE_URL = env("SITE_URL")
CSRF_TRUSTED_ORIGINS = [SITE_URL]
ROOT_URLCONF = "backend.application.urls"
WSGI_APPLICATION = "cellarium.nexus.backend.application.wsgi.application"
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"
ENVIRONMENT = env("ENVIRONMENT", default="local")

# Application Definition
# ---------------------
INSTALLED_APPS = [
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "unfold",
    "unfold.contrib.filters",
    "django.contrib.admin",
    "rest_framework",
    "storages",
    "django_json_widget",
    "cellarium.nexus.backend.core",
    "cellarium.nexus.backend.cell_management",
    "cellarium.nexus.backend.curriculum",
    "cellarium.nexus.backend.ingest_management",
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

# Template Configuration
# ----------------------
TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [BASE_DIR / "cellarium/nexus/backend/templates"],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

default_db_config = {
    "ENGINE": "django.db.backends.postgresql",
    "NAME": env("DB_NAME"),
    "USER": env("DB_USER"),
    "PASSWORD": env("DB_PASSWORD"),
}

if ENVIRONMENT != "local":
    # Cloud Run environment using Cloud SQL socket
    default_db_config["HOST"] = f'/cloudsql/{env("DB_INSTANCE_CONNECTION_NAME")}'
else:
    # Local or other environments using standard host/port
    default_db_config["HOST"] = env("DB_HOST", default="localhost")  # Add default for local
    default_db_config["PORT"] = env("DB_PORT", default="5432")  # Add default for local

DATABASES = {"default": default_db_config}

# Authentication and Security
# --------------------------
AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]

# Internationalization and Localization
# ----------------------------------
LANGUAGE_CODE = "en-us"
TIME_ZONE = "UTC"
USE_I18N = True
USE_TZ = True

# Django Unfold Admin Settings
# --------------------------
UNFOLD = {
    "SITE_TITLE": "Cellarium Nexus",
    "SITE_HEADER": "Cellarium Nexus Admin",
    "SITE_SUBHEADER": "Data Management Interface",
    # Logo configuration
    "SITE_LOGO": {
        "light": lambda request: static("cellarium-nexus-logo-light.svg"),
        "dark": lambda request: static("cellarium-nexus-logo-dark.svg"),
    },
    # Favicon configuration
    "SITE_FAVICONS": [
        {
            "rel": "icon",
            "sizes": "32x32",
            "href": lambda request: static("favicon.svg"),
            "type": "image/svg+xml",
        },
    ],
    # UI preferences
    "SHOW_HISTORY": True,
    "SHOW_VIEW_ON_SITE": True,
    "SHOW_BACK_BUTTON": True,
    "COLORS": {
        "base": {
            "50": "250 250 252",
            "100": "244 244 247",
            "200": "228 228 236",
            "300": "208 208 220",
            "400": "150 150 170",
            "500": "107 107 128",
            "600": "70 70 95",
            "700": "50 50 75",
            "800": "30 30 50",
            "900": "15 15 30",
            "950": "8 8 20",
        },
        "primary": {
            "50": "240 238 255",
            "100": "225 220 255",
            "200": "190 180 255",
            "300": "150 130 255",
            "400": "120 90 255",
            "500": "90 50 255",
            "600": "70 30 220",
            "700": "55 20 180",
            "800": "40 15 140",
            "900": "30 10 100",
            "950": "20 5 70",
        },
        "font": {
            "subtle-light": "var(--color-base-500)",
            "subtle-dark": "var(--color-base-400)",
            "default-light": "var(--color-base-700)",
            "default-dark": "var(--color-base-300)",
            "important-light": "var(--color-base-900)",
            "important-dark": "var(--color-base-100)",
        },
    },
    # Sidebar configuration
    "SIDEBAR": {
        "show_search": True,
        "show_all_applications": True,
        "navigation": [
            {
                "title": _("Omics Data"),
                "separator": True,
                "collapsible": False,
                "items": [
                    {
                        "title": _("Cell Information"),
                        "icon": "science",
                        "link": reverse_lazy("admin:cell_management_cellinfo_changelist"),
                    },
                    {
                        "title": _("Curriculums"),
                        "icon": "school",
                        "link": reverse_lazy("admin:curriculum_curriculum_changelist"),
                    },
                    {
                        "title": _("Feature Schemas"),
                        "icon": "schema",
                        "link": reverse_lazy("admin:cell_management_featureschema_changelist"),
                    },
                ],
            },
            {
                "title": _("Ingest Management"),
                "separator": True,
                "collapsible": True,
                "expanded": False,
                "items": [
                    {
                        "title": _("BigQuery Datasets"),
                        "icon": "database",
                        "link": reverse_lazy("admin:cell_management_bigquerydataset_changelist"),
                    },
                    {
                        "title": _("Data Ingests"),
                        "icon": "upload_file",
                        "link": reverse_lazy("admin:ingest_management_ingestinfo_changelist"),
                    },
                    {
                        "title": _("Column Mappers"),
                        "icon": "table_rows",
                        "link": reverse_lazy("admin:ingest_management_columnmapping_changelist"),
                    },
                    {
                        "title": _("Validation Reports"),
                        "icon": "fact_check",
                        "link": reverse_lazy("admin:ingest_management_validationreport_changelist"),
                    },
                ],
            },
            {
                "title": _("Administration"),
                "separator": True,
                "collapsible": True,
                "items": [
                    {
                        "title": _("Users"),
                        "icon": "people",
                        "link": reverse_lazy("admin:auth_user_changelist"),
                        "permission": lambda request: request.user.is_superuser,
                    },
                    {
                        "title": _("Groups"),
                        "icon": "groups",
                        "link": reverse_lazy("admin:auth_group_changelist"),
                        "permission": lambda request: request.user.is_superuser,
                    },
                ],
            },
        ],
    },
    "TABS": [
        {
            "models": [
                "cell_management.cellinfo",
                "cell_management.cellfeatureinfo",
            ],
            "items": [
                {
                    "title": _("Cells"),
                    "link": reverse_lazy("admin:cell_management_cellinfo_changelist"),
                },
                {
                    "title": _("Cell Features"),
                    "link": reverse_lazy("admin:cell_management_cellfeatureinfo_changelist"),
                },
            ],
        },
    ],
}

# Google Cloud Storage Settings
# --------------------------
GCP_PROJECT_ID = env("GCP_PROJECT_ID")
GCP_APPLICATION_BILLING_LABEL = env("GCP_APPLICATION_BILLING_LABEL")
BUCKET_NAME_PRIVATE = env("BUCKET_NAME_PRIVATE")
BUCKET_NAME_PUBLIC = env("BUCKET_NAME_PUBLIC")
STATIC_LOCATION = "static-files"

# Static files settings
# --------------------
STATIC_URL = f"https://storage.googleapis.com/{BUCKET_NAME_PUBLIC}/{STATIC_LOCATION}/"
STATICFILES_DIRS = [BASE_DIR / "cellarium/nexus/backend/application/static"]

# Google Cloud Storage settings
GS_QUERYSTRING_AUTH = False  # Don't use signed URLs for static files
GS_BUCKET_NAME = BUCKET_NAME_PUBLIC
GS_PROJECT_ID = GCP_PROJECT_ID

# Don't set ACLs (using uniform bucket-level access)
GS_DEFAULT_ACL = None
GS_OBJECT_PARAMETERS = {"cache-control": "public, max-age=86400"}  # Cache for 24 hours

# CORS configuration for GCS
GS_CORS_ORIGIN_ALLOW_ALL = True
GS_CORS_METHODS = ["GET", "OPTIONS"]
GS_CORS_MAX_AGE = 3600  # 1 hour cache
GS_CORS_HEADERS = [
    "Content-Type",
    "Access-Control-Allow-Origin",
    "Access-Control-Allow-Headers",
    "Origin",
    "X-Requested-With",
]

STORAGES = {
    "default": {
        "BACKEND": "storages.backends.gcloud.GoogleCloudStorage",
        "OPTIONS": {
            "project_id": GCP_PROJECT_ID,
            "bucket_name": BUCKET_NAME_PRIVATE,
            "location": "",
            "querystring_auth": False,
        },
    },
    "staticfiles": {
        "BACKEND": "storages.backends.gcloud.GoogleCloudStorage",
        "OPTIONS": {
            "project_id": GCP_PROJECT_ID,
            "bucket_name": BUCKET_NAME_PUBLIC,
            "location": STATIC_LOCATION,
            "querystring_auth": False,
        },
    },
}

BACKEND_PIPELINE_DIR = "pipeline"
INGEST_INPUT_FILE_MAX_SIZE = 3e9  # 3 Gb
MAX_ADATA_FILES_PER_VALIDATION_BATCH = 50  # Maximum number of AnnData files per validation batch

# Vertex AI Pipeline settings
PIPELINE_SERVICE_ACCOUNT = env("PIPELINE_SERVICE_ACCOUNT", default=None)
PIPELINE_BASE_IMAGE = env("PIPELINE_BASE_IMAGE")
PIPELINE_ROOT_PATH = env("PIPELINE_ROOT_PATH", default=None)
