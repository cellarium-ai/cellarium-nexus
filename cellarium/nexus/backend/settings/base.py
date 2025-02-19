from django.templatetags.static import static
from django.urls import reverse_lazy
from django.utils.translation import gettext_lazy as _
from nexus.backend.settings import BASE_DIR, env

SECRET_KEY = env("SECRET_KEY")

DEBUG = True

ALLOWED_HOSTS = []
SITE_URL = "http://localhost:8000"

# Application definition

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
    "nexus.backend.app_modules.cell_management",
    "nexus.backend.app_modules.core",
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

ROOT_URLCONF = "backend.urls"

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

WSGI_APPLICATION = "nexus.backend.wsgi.application"

# Database
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": env("DB_NAME"),
        "USER": env("DB_USER"),
        "PASSWORD": env("DB_PASSWORD"),
        "HOST": env("DB_HOST"),
        "PORT": env("DB_PORT"),
    }
}
# Password validation
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

# Internationalization
LANGUAGE_CODE = "en-us"

TIME_ZONE = "UTC"

USE_I18N = True

USE_TZ = True


# Default primary key field type
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

# Django Unfold admin settings
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
                "title": _("Data Management"),
                "separator": True,
                "collapsible": False,
                "items": [
                    {
                        "title": _("Cell Information"),
                        "icon": "science",
                        "link": reverse_lazy("admin:cell_management_cellinfo_changelist"),
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
                        "link": reverse_lazy("admin:cell_management_ingestinfo_changelist"),
                    },
                    {
                        "title": _("Column Mappers"),
                        "icon": "table_rows",
                        "link": reverse_lazy("admin:cell_management_columnmapping_changelist"),
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

# Static files settings
# Static files (CSS, JavaScript, Images)
STATIC_URL = "static/"
STATIC_ROOT = BASE_DIR / "staticfiles"
STATICFILES_DIRS = [
    BASE_DIR / "cellarium/nexus/backend/static",
]

# GS_BUCKET_NAME = "cellarium-file-system"

# STORAGES = {"default": {"BACKEND": "storages.backends.gcloud.GoogleCloudStorage", "OPTIONS": {}}}
GCP_PROJECT_ID = "dsp-cell-annotation-service"
BUCKET_NAME_PRIVATE = "cellarium-nx-file-system"
BUCKET_NAME_PUBLIC = "cellarium-nx-public-files"
STATIC_LOCATION = "static-files"
# STORAGES = {
#     "default": {
#         "BACKEND": "nexus.backend.app_storages_backend.CustomGoogleCloudStorage",
#         "OPTIONS": {
#             "bucket_name": BUCKET_NAME_PRIVATE,
#             "default_acl": "publicRead",
#             "file_overwrite": False,
#         },
#     },
#     "staticfiles": {
#         "BACKEND": "nexus.backend.app_storages_backend.CustomGoogleCloudStorage",
#         "OPTIONS": {
#             "bucket_name": BUCKET_NAME_PUBLIC,
#             "location": STATIC_LOCATION,
#             "default_acl": "publicRead",
#             "file_overwrite": True,
#         },
#     },
# }
# STATIC_URL = f"/static/"
# STATIC_URL = f"https://storage.googleapis.com/{BUCKET_NAME_PUBLIC}/{STATIC_LOCATION}/"
BACKEND_PIPELINE_DIR = "pipeline"
# STATIC_ROOT = "django-static"
# DEFAULT_FILE_STORAGE = "storages.backends.gcloud.GoogleCloudStorage"
# STATICFILES_DIRS = [
#     BASE_DIR / "cellarium/nexus/assets",  # Include the assets directory
# ]
# STATIC_URL = f"https://storage.googleapis.com/{GS_BUCKET_NAME}/static/"
#
# # Optional: Media files (if separate from static files)
# MEDIA_URL = f"https://storage.googleapis.com/{GS_BUCKET_NAME}/media/"
# MEDIA_ROOT = "media/"
#
# # Static files storage (optional, if different from media files)
# STATICFILES_STORAGE = "storages.backends.gcloud.GoogleCloudStorage"
