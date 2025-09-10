"""
Constants used across the admin interface.
"""

from django.utils.translation import gettext_lazy as _

# BigQuery related messages and formats
BIGQUERY_SUCCESS_MESSAGE_LINK_FORMAT = (
    '<a href="{url_link}" target="_blank" style="text-decoration: underline;">View in BigQuery</a>'
)
BIGQUERY_SUCCESS_MESSAGE_TEXT = _("BigQuery dataset in GCP was created successfully.")

# Template paths
CHANGELIST_ACTION_FORM = "admin/custom_templates/changelist_action_with_form.html"

# CSV related constants
REQUIRED_CSV_FILE_COLUMNS = ["gcs_file_path"]

# Form titles and messages
PREPARE_EXTRACT_TABLES_TITLE = _("Extract Curriculum Job")
PREPARE_BUTTON_TITLE = _("Submit")
NO_DATASETS_ERROR = _("No BigQuery datasets available. You need to create one to be able to extract data.")
MULTIPLE_DATASETS_ERROR = _("Multiple BigQuery datasets exist. Please select a BigQuery dataset for extraction.")

# Pipeline success messages
EXTRACT_PIPELINE_SUCCESS_MESSAGE = _(
    'Extract curriculum job started successfully. <a href="{0}" target="_blank" style="text-decoration: underline;">Vertex AI Dashboard</a>'
)

# Workflows related constants
BINS_PER_WORKER = 32

# CellInfo schema columns for extraction
CELL_INFO_EXTRACT_COLUMNS_CHOICES = [
    ("original_id", "Original Id"),
    ("donor_id", "Donor Id"),
    ("cell_type", "Cell Type"),
    ("assay", "Assay"),
    ("development_stage", "Development Stage"),
    ("tissue", "Tissue"),
    ("disease", "Disease"),
    ("organism", "Organism"),
    ("self_reported_ethnicity", "Self Reported Ethnicity"),
    ("sex", "Sex"),
    ("suspension_type", "Suspension Type"),
    ("total_mrna_umis", "Total Mrna Umis"),
    ("cell_type_ontology_term_id", "Cell Type Ontology Term Id"),
    ("assay_ontology_term_id", "Assay Ontology Term Id"),
    ("development_stage_ontology_term_id", "Development Stage Ontology Term Id"),
    ("tissue_ontology_term_id", "Tissue Ontology Term Id"),
    ("disease_ontology_term_id", "Disease Ontology Term Id"),
    ("organism_ontology_term_id", "Organism Ontology Term Id"),
    ("self_reported_ethnicity_ontology_term_id", "Self Reported Ethnicity Ontology Term Id"),
    ("sex_ontology_term_id", "Sex Ontology Term Id"),
    ("tag", "Tag"),
]

CELL_INFO_EXTRACT_BIN_KEYS_COLUMNS_CHOICES = [
    ("donor_id", "Donor Id"),
    ("cell_type", "Cell Type"),
    ("assay", "Assay"),
    ("development_stage", "Development Stage"),
    ("tissue", "Tissue"),
    ("disease", "Disease"),
    ("organism", "Organism"),
    ("self_reported_ethnicity", "Self Reported Ethnicity"),
    ("sex", "Sex"),
    ("suspension_type", "Suspension Type"),
    ("cell_type_ontology_term_id", "Cell Type Ontology Term Id"),
    ("assay_ontology_term_id", "Assay Ontology Term Id"),
    ("development_stage_ontology_term_id", "Development Stage Ontology Term Id"),
    ("tissue_ontology_term_id", "Tissue Ontology Term Id"),
    ("disease_ontology_term_id", "Disease Ontology Term Id"),
    ("organism_ontology_term_id", "Organism Ontology Term Id"),
    ("self_reported_ethnicity_ontology_term_id", "Self Reported Ethnicity Ontology Term Id"),
    ("sex_ontology_term_id", "Sex Ontology Term Id"),
    ("tag", "Tag"),
]
