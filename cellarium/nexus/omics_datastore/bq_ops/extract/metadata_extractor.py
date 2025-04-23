"""
Extract metadata from BigQuery tables for CAS data extracts.
"""

import json
import logging
import os
import ssl
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import certifi
import numpy as np
import pandas as pd
from google.cloud import bigquery
from smart_open import open

from cellarium.nexus.omics_datastore import bq_sql

# Configure SSL context for aiohttp
ssl_context = ssl.create_default_context(cafile=certifi.where())
os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()
os.environ["SSL_CERT_FILE"] = certifi.where()


logger = logging.getLogger(__name__)

# Template paths
TEMPLATE_DIR = Path(__file__).parent.parent.parent / "sql_templates" / "metadata"
EXTRACT_METADATA_TEMPLATE = TEMPLATE_DIR / "extract_metadata.sql.mako"
MEASURED_GENES_TEMPLATE = TEMPLATE_DIR / "measured_genes.sql.mako"

COLUMNS_TO_IGNORE = ["id", "extract_bin", "original_id", "metadata_extra"]


@dataclass
class ExtractMetadata:
    """
    Store metadata about the extract.

    :param total_bins: Total number of extract bins
    :param last_bin_size: Size of the last bin
    :param filters: Filters used for the extract
    """

    total_bins: int
    last_bin_size: int
    filters: dict[str, Any] | None


class MetadataExtractor:
    """
    Extract and save metadata from BigQuery tables.

    :param client: BigQuery client
    :param project: GCP project ID
    :param dataset: BigQuery dataset ID
    :param extract_table_prefix: Prefix for extract tables
    :param filters: Query filters
    """

    def __init__(
        self,
        client: bigquery.Client,
        project: str,
        dataset: str,
        extract_table_prefix: str,
        filters: dict[str, Any] | None = None,
    ) -> None:
        self.client = client
        self.project = project
        self.dataset = dataset
        self.prefix = extract_table_prefix
        self.filters = filters

    def execute_query(self, sql: str) -> Any:
        """
        Execute BigQuery SQL query.

        :param sql: SQL query string to execute

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: Query results
        """
        query = self.client.query(sql)
        return query.result()

    def get_extract_metadata(self) -> ExtractMetadata:
        """
        Get metadata about the extract bins.

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: ExtractMetadata object with bin information
        """
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
        )
        sql = bq_sql.render(str(EXTRACT_METADATA_TEMPLATE), template_data)
        result = list(self.execute_query(sql))[0]

        return ExtractMetadata(total_bins=result.total_bins, last_bin_size=result.last_bin_size, filters=self.filters)

    def get_categorical_columns(self) -> dict[str, list[str]]:
        """
        Get categorical column metadata.

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: Dictionary mapping column names to their unique values
        """
        table_ref = f"{self.project}.{self.dataset}.{self.prefix}__extract_cell_info"
        table = self.client.get_table(table_ref)

        categorical_columns = {}
        for field in table.schema:
            if field.field_type == "STRING" and field.name not in COLUMNS_TO_IGNORE:
                template_data = bq_sql.TemplateData(
                    project=self.project, dataset=self.dataset, extract_table_prefix=self.prefix, column_name=field.name
                )
                sql = bq_sql.render(str(TEMPLATE_DIR / "categorical_values.sql.mako"), template_data)
                results = self.execute_query(sql)
                categorical_columns[field.name] = [row[0] for row in results]

        return categorical_columns

    def get_measured_genes_info(self) -> pd.DataFrame:
        """
        Get information about measured genes grouped by tag.

        Creates a matrix where rows are tags and columns are features, with values
        indicating whether a feature is measured for that tag.

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: DataFrame with gene measurement matrix (tags x features)
        """
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
        )
        sql = bq_sql.render(str(MEASURED_GENES_TEMPLATE), template_data)
        query_result = self.client.query(sql).result()

        # Convert query results to DataFrame
        all_features = set()
        gene_measurement_data = []
        tags = []
        feature_info = {}  # Store feature info for labels

        for row in query_result:
            tags.append(row.tag)
            # Create a dictionary mapping feature IDs to True
            measured_features = {}
            for feature in row.features:
                # Access dictionary fields using string keys
                feature_id = feature["id"]
                measured_features[feature_id] = True
                # Store feature info for labels if we haven't seen it yet
                if feature_id not in feature_info:
                    feature_info[feature_id] = (feature["symbol"], feature["ensemble_id"])
            # Keep track of all unique feature IDs
            all_features.update(measured_features.keys())
            gene_measurement_data.append(measured_features)

        # Create measurement matrix
        all_features = sorted(all_features)
        gene_expression_mask = np.zeros((len(tags), len(all_features)))

        for i, measured_features in enumerate(gene_measurement_data):
            for j, feature_id in enumerate(all_features):
                gene_expression_mask[i, j] = measured_features.get(feature_id, False)

        # Create feature labels using symbol and ensemble_id
        feature_labels = [f"{feature_info[feature_id][0]}_{feature_info[feature_id][1]}" for feature_id in all_features]

        return pd.DataFrame(data=gene_expression_mask, index=pd.Index(tags, name="tag"), columns=feature_labels)

    def save_metadata(
        self,
        bucket_name: str,
        extract_bucket_path: str,
        metadata_dir: str = "shared_metadata",
    ) -> None:
        """
        Save all extract metadata to a single JSON file in GCS bucket.

        :param bucket_name: GCS bucket name
        :param extract_bucket_path: Path within bucket for extract
        :param metadata_dir: Directory name for metadata files

        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail
        :raise IOError: If file operations fail
        """
        base_path = f"gs://{bucket_name}/{extract_bucket_path}/{metadata_dir}"

        # Collect all metadata
        logger.info("Extracting extract metadata...")
        extract_meta = self.get_extract_metadata()

        logger.info("Extracting categorical columns metadata...")
        categorical_meta = self.get_categorical_columns()

        logger.info("Extracting measured genes info...")
        genes_info = self.get_measured_genes_info()

        # Combine all metadata into a single dictionary
        combined_metadata = {
            "general_info": extract_meta.__dict__,
            "category_metadata": categorical_meta,
            "measured_genes_mask": genes_info.to_dict(orient="records"),
        }

        # Save combined metadata to a single file
        logger.info("Saving combined metadata...")
        with open(f"{base_path}/extract_metadata.json", "w") as f:
            json.dump(combined_metadata, f, indent=2)
