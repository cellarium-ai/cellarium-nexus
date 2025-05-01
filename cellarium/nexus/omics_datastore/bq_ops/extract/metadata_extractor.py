"""
Extract metadata from BigQuery tables for CAS data extracts.
"""

import json
import logging
import os
import ssl
from pathlib import Path
from typing import Any

import certifi
import numpy as np
import pandas as pd
from google.cloud import bigquery
from smart_open import open

from cellarium.nexus.omics_datastore import bq_sql
from cellarium.nexus.shared.schemas.omics_datastore import ExtractMetadata

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
        extract_bin_size: int | None = None,
    ) -> None:
        self.client = client
        self.project = project
        self.dataset = dataset
        self.prefix = extract_table_prefix
        self.filters = filters
        self.extract_bin_size = extract_bin_size

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

        return ExtractMetadata(
            total_bins=result.total_bins,
            last_bin_size=result.last_bin_size,
            total_cells=result.total_cells,
            filters=self.filters,
            extract_bin_size=self.extract_bin_size,
        )

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

    def compose_extract_metadata(self) -> ExtractMetadata:
        """
        Create and return a complete ExtractMetadata object with all metadata.

        :raise google.api_core.exceptions.GoogleAPIError: If BigQuery operations fail

        :return: The complete ExtractMetadata object
        """
        logger.info("Extracting extract metadata...")
        extract_meta = self.get_extract_metadata()

        logger.info("Extracting categorical columns metadata...")
        categorical_meta = self.get_categorical_columns()

        # Update the extract_meta with additional attributes
        extract_meta.category_metadata = categorical_meta

        return extract_meta
