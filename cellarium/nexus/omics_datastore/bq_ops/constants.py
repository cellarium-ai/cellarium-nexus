BQ_INGEST_TABLE_NAME = "ingest_info"
BQ_FEATURE_INFO_TABLE_NAME = "feature_info"
BQ_CELL_INFO_TABLE_NAME = "cell_info"
BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME = "raw_matrix_coo"

# Extract table names
BQ_EXTRACT_FEATURE_INFO_TABLE_NAME = "__extract_feature_info"
BQ_EXTRACT_CELL_INFO_TABLE_NAME = "__extract_cell_info"
BQ_EXTRACT_CELL_INFO_RANDOMIZED_TABLE_NAME = "__extract_cell_info_randomized"
BQ_EXTRACT_MATRIX_COO_TABLE_NAME = "__extract_matrix_coo"

INGEST_INGEST_FILE_NAME = "ingest-info.avro"
INGEST_CELL_INFO_FILE_NAME = "cell-info.avro"
INGEST_FEATURE_INFO_FILE_NAME = "feature-info.avro"
INGEST_RAW_COUNTS_FILE_NAME_FORMAT = "raw-counts-{batch_index:06}.csv"
INGEST_RAW_COUNTS_FILE_PATTERN = "raw-counts-*.csv"


# Column Names
OBS_CELL_INFO_ORIGINAL_ID = "original_id"
OBS_NEXUS_ID = "id"
OBS_INGEST_ID = "ingest_id"
OBS_TOTAL_MRNA_UMIS = "total_mrna_umis"
OBS_METADATA_EXTRA = "metadata_extra"
OBS_BIGQUERY_DATASET = "bigquery_dataset"
OBS_TAG = "tag"
VAR_FEATURE_INFO_ORIGINAL_ID = "original_id"
VAR_NEXUS_ID = "id"
VAR_INGEST_ID = "ingest_id"
VAR_METADATA_EXTRA = "metadata_extra"
VAR_TAG = "tag"

# Metadata storage paths
SHARED_META_DIR_NAME = "shared_metadata"
MEASURED_GENES_INFO_FILE_NAME = "measured_genes_info.csv"
CATEGORICAL_COLUMNS_META_FILE_NAME = "categorical_columns_metadata.pkl"
