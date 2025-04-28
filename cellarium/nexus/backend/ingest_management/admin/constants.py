from django.utils.translation import gettext_lazy as _

BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT = (
    '<a href="{url_link}" target="_blank" style="text-decoration: underline;">View in BigQuery</a>'
)
BIGQUERY_SUCCESS_MESSAGE_TEXT = _("BigQuery dataset in GCP was created successfully.")
CHANGELIST_ACTION_FORM = "admin/custom_templates/changelist_action_with_form.html"
GCS_PATH_COLUMN = "gcs_file_path"
TAGS_COLUMN = "tag"
REQUIRED_CSV_FILE_COLUMNS = [GCS_PATH_COLUMN, TAGS_COLUMN]

# Form titles and messages
PREPARE_EXTRACT_TABLES_TITLE = _("Prepare Extract Tables")
PREPARE_BUTTON_TITLE = _("Prepare")
NO_DATASETS_ERROR = _("No BigQuery datasets available")
MULTIPLE_DATASETS_ERROR = _("Multiple BigQuery datasets exist. Please select a BigQuery dataset for extraction.")
EXTRACT_SUCCESS_MESSAGE = _("Extract tables prepared successfully")

# Pipeline success messages
INGEST_PIPELINE_SUCCESS_MESSAGE = _(
    'Data ingestion pipeline started successfully. <a href="{0}" target="_blank" style="text-decoration: underline;">Vertex AI Dashboard</a>'
)
