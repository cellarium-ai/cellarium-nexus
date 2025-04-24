"""
Extract data from BigQuery tables into AnnData files.
"""

import concurrent.futures as concurrency
import logging
import multiprocessing
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from google.cloud import bigquery
from google.cloud.bigquery_storage import BigQueryReadClient, types
from scipy.sparse import coo_matrix
from tenacity import before_log, retry, stop_after_attempt, wait_exponential

from cellarium.nexus.omics_datastore import bq_sql
from cellarium.nexus.omics_datastore.bq_ops import constants
from cellarium.nexus.shared import schemas

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Template paths
TEMPLATE_DIR = Path(__file__).parent.parent.parent / "sql_templates" / "extract"
GET_CELLS_IN_BIN_RANGE_TEMPLATE = TEMPLATE_DIR / "get_cells_in_bin_range.sql.mako"
GET_FEATURES_TEMPLATE = TEMPLATE_DIR / "get_features.sql.mako"


class DataExtractor:
    """
    Extract data from BigQuery tables into AnnData files.
    """

    def __init__(
        self,
        *,
        client: bigquery.Client,
        project: str,
        dataset: str,
        extract_table_prefix: str,
    ) -> None:
        """
        Initialize data extractor.

        :param client: BigQuery client instance
        :param project: GCP project ID
        :param dataset: BigQuery dataset ID
        :param extract_table_prefix: Prefix for extract table names
        """
        self.client = client
        self.project = project
        self.dataset = dataset
        self.prefix = extract_table_prefix
        self.read_client = BigQueryReadClient()

    def execute_query(self, *, sql: str) -> Any:
        """
        Execute a BigQuery SQL query.

        :param sql: SQL query string to execute

        :raise google.api_core.exceptions.GoogleAPIError: If query execution fails

        :return: Query results
        """
        query = self.client.query(sql)
        return query.result()

    def get_features(self) -> list[schemas.FeatureSchema]:
        """
        Get features from extract feature info table.

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: List of schemas.FeatureSchema objects
        """
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
        )
        sql = bq_sql.render(str(GET_FEATURES_TEMPLATE), template_data)
        result = self.execute_query(sql=sql)
        return [
            schemas.FeatureSchema(id=row["id"], symbol=row["symbol"], ensemble_id=row["ensemble_id"]) for row in result
        ]

    def get_cells_in_bin_range(
        self,
        *,
        start_bin: int,
        end_bin: int,
        obs_columns: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        """
        Get cells within specified bin range.

        :param start_bin: Starting bin number (inclusive)
        :param end_bin: Ending bin number (inclusive)
        :param obs_columns: Optional observation columns to include

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: Sequence of cell data dictionaries
        """
        template_data = bq_sql.TemplateData(
            project=self.project,
            dataset=self.dataset,
            extract_table_prefix=self.prefix,
            select=obs_columns,
            start_bin=start_bin,
            end_bin=end_bin,
        )
        sql = bq_sql.render(str(GET_CELLS_IN_BIN_RANGE_TEMPLATE), template_data)
        return [x for x in self.execute_query(sql=sql)]

    def get_matrix_data(
        self,
        *,
        start_bin: int,
        end_bin: int,
    ) -> Any:
        """
        Get matrix data for specified bin range using storage API.

        :param start_bin: Starting bin number (inclusive)
        :param end_bin: Ending bin number (inclusive)

        :raise google.api_core.exceptions.GoogleAPIError: If query fails

        :return: Matrix data rows
        """
        table = f"projects/{self.project}/datasets/{self.dataset}/tables/{self.prefix}{constants.BQ_EXTRACT_MATRIX_COO_TABLE_NAME}"

        requested_session = types.ReadSession()
        requested_session.table = table
        requested_session.data_format = types.DataFormat.AVRO

        requested_session.read_options.selected_fields = ["cell_id", "feature_data"]
        requested_session.read_options.row_restriction = f"extract_bin BETWEEN {start_bin} AND {end_bin}"

        parent = f"projects/{self.project}"
        session = self.read_client.create_read_session(
            parent=parent,
            read_session=requested_session,
            max_stream_count=1,
        )

        logger.info(f"Estimated bytes to scan: {session.estimated_total_bytes_scanned}")
        reader = self.read_client.read_rows(session.streams[0].name)
        return reader.rows(session)

    def convert_matrix_data_to_coo(
        self,
        *,
        matrix_data: Any,
        cell_index_to_row_num: dict[int, int],
        feature_id_to_col_num: dict[int, int],
    ) -> tuple[list[int], list[int], list[float]]:
        """
        Convert matrix data to COO format components.

        :param matrix_data: Raw matrix data from storage API
        :param cell_index_to_row_num: Mapping of cell indices to row numbers
        :param feature_id_to_col_num: Mapping of feature indices to column numbers

        :raise KeyError: If index mappings are invalid

        :return: Tuple of (rows, columns, data) for COO matrix
        """
        rows = []
        columns = []
        data = []

        for row in matrix_data:
            cell_id = row["cell_id"]
            feature_data = row["feature_data"]

            try:
                row_num = cell_index_to_row_num[cell_id]
            except KeyError as exc:
                raise KeyError(f"Invalid cell id: {cell_id}") from exc

            for entry in feature_data:
                feature_index = entry["feature_id"]
                raw_counts = entry["raw_counts"]

                try:
                    col_num = feature_id_to_col_num[feature_index]
                except KeyError as exc:
                    raise KeyError(f"Invalid feature_index: {feature_index}") from exc

                rows.append(row_num)
                columns.append(col_num)
                data.append(raw_counts)

        return rows, columns, data

    def extract_bin_to_anndata(
        self,
        *,
        bin_number: int,
        output_path: Path,
        extract_metadata: schemas.ExtractMetadata,
        obs_columns: list[str] | None = None,
    ) -> None:
        """
        Extract a single bin to an AnnData file.

        :param bin_number: Bin number to extract
        :param output_path: Local path to save AnnData file
        :param extract_metadata: ExtractMetadata instance with metadata for the extract
        :param obs_columns: Optional observation columns to include

        :raise google.api_core.exceptions.GoogleAPIError: If queries fail
        :raise IOError: If file operations fail
        """
        logger.info(f"Extracting bin {bin_number}...")

        # Get features
        features = self.get_features()
        feature_ids = []
        feature_id_to_col_num = {}

        for col_num, feature in enumerate(features):
            feature_ids.append(feature.ensemble_id)  # Use ensemble_id as feature identifier
            feature_id_to_col_num[feature.id] = col_num

        # Get cells
        cells = self.get_cells_in_bin_range(start_bin=bin_number, end_bin=bin_number, obs_columns=obs_columns)
        original_cell_ids = []
        ingest_ids = []
        obs_columns_values = {k: [] for k in (obs_columns or [])}
        cell_id_to_row_num = {}

        for row_num, cell in enumerate(cells):
            original_cell_ids.append(cell["id"])
            ingest_ids.append(cell["ingest_id"])
            cell_id_to_row_num[cell["id"]] = row_num
            for obs_column in obs_columns or []:
                obs_columns_values[obs_column].append(cell[obs_column])

        # Get matrix data
        matrix_data = self.get_matrix_data(start_bin=bin_number, end_bin=bin_number)

        rows, columns, data = self.convert_matrix_data_to_coo(
            matrix_data=matrix_data,
            cell_index_to_row_num=cell_id_to_row_num,
            feature_id_to_col_num=feature_id_to_col_num,
        )

        # Create AnnData object
        counts = coo_matrix((data, (rows, columns)), shape=(len(cells), len(features)), dtype=np.float32)
        adata = ad.AnnData(counts.tocsr())
        adata.obs.index = pd.Index(original_cell_ids)
        adata.var.index = pd.Index(feature_ids)

        # Add observation columns
        for obs_column, values in obs_columns_values.items():
            # Convert to pandas Series first
            series = pd.Series(values, index=adata.obs.index)

            # Handle None values based on data type
            if series.dropna().empty:
                # Skip columns with all None values
                logger.info(f"Skipping column '{obs_column}' as all values are None")
                continue

            # For numeric-like columns, None becomes NaN
            # For string-like columns, None becomes empty string
            if pd.api.types.is_numeric_dtype(series.dropna()):
                # Replace None with NaN for numeric columns
                adata.obs[obs_column] = series
            else:
                # Replace None with empty string for string columns
                adata.obs[obs_column] = series.fillna("")

        # Use the categorical columns from the provided metadata
        categorical_columns = extract_metadata.category_metadata
        logger.info(f"Using provided metadata with {len(categorical_columns)} categorical columns")

        # Apply categorical columns to the AnnData object
        for column, categories in categorical_columns.items():
            if column in adata.obs:
                adata.obs[column] = pd.Categorical(values=pd.Series(adata.obs[column]), categories=categories)

        # Save file
        adata.write_h5ad(output_path, compression="gzip")
        logger.info(f"Saved bin {bin_number} to {output_path}")


def child_init(log_level: str) -> None:
    """
    Configure logging for child processes.

    :param log_level: Logging level to set
    """
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] [%(processName)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    before=before_log(logger, logging.INFO),
)
def perform_extraction(
    *,
    extractor: DataExtractor,
    bin_number: int,
    output_path: Path,
    extract_metadata: schemas.ExtractMetadata,
    obs_columns: list[str] | None = None,
) -> None:
    """
    Perform extraction with retry logic.

    :param extractor: DataExtractor instance
    :param bin_number: Bin number to extract
    :param output_path: Local path to save AnnData file
    :param obs_columns: Optional observation columns to include
    :param extract_metadata: Optional pre-loaded ExtractMetadata instance

    :raise: Will raise any exception after retry attempts are exhausted
    """
    extractor.extract_bin_to_anndata(
        bin_number=bin_number,
        output_path=output_path,
        extract_metadata=extract_metadata,
        obs_columns=obs_columns,
    )


def extract_bin_to_anndata_worker(
    *,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    bin_number: int,
    output_path: Path,
    extract_metadata: schemas.ExtractMetadata,
    obs_columns: list[str] | None = None,
) -> None:
    """
    Worker function to extract a single bin to an AnnData file.
    Creates its own BigQuery clients to avoid pickling issues.

    :param project: GCP project ID
    :param dataset: BigQuery dataset ID
    :param extract_table_prefix: Prefix for extract table names
    :param bin_number: Bin number to extract
    :param output_path: Local path to save AnnData file
    :param extract_metadata: ExtractMetadata instance with metadata for the extract
    :param obs_columns: Optional observation columns to include

    :raise google.api_core.exceptions.GoogleAPIError: If queries fail
    :raise IOError: If file operations fail
    """
    # Create clients in the worker process
    client = bigquery.Client(project=project)
    extractor = DataExtractor(
        client=client,
        project=project,
        dataset=dataset,
        extract_table_prefix=extract_table_prefix,
    )

    try:
        perform_extraction(
            extractor=extractor,
            bin_number=bin_number,
            output_path=output_path,
            extract_metadata=extract_metadata,
            obs_columns=obs_columns,
        )
    finally:
        client.close()


def extract_bins(
    *,
    client: bigquery.Client,
    project: str,
    dataset: str,
    extract_table_prefix: str,
    bins: list[int],
    output_dir: Path,
    extract_metadata: schemas.ExtractMetadata,
    obs_columns: list[str] | None = None,
    max_workers: int | None = None,
) -> None:
    """
    Extract multiple bins in parallel.

    :param client: BigQuery client instance (used only for validation)
    :param project: GCP project ID
    :param dataset: BigQuery dataset ID
    :param extract_table_prefix: Prefix for extract table names
    :param bins: List of bin numbers to extract
    :param output_dir: Local directory to save AnnData files
    :param extract_metadata: ExtractMetadata instance with metadata for the extract
    :param obs_columns: Optional observation columns to include
    :param max_workers: Maximum number of parallel workers

    :raise google.api_core.exceptions.GoogleAPIError: If queries fail
    :raise IOError: If file operations fail
    """
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    with concurrency.ProcessPoolExecutor(
        max_workers=max_workers, initializer=child_init, initargs=(logging.getLevelName(logging.INFO),)
    ) as executor:
        futures = []
        for bin_num in bins:
            output_path = output_dir / f"extract_{bin_num}.h5ad"
            future = executor.submit(
                extract_bin_to_anndata_worker,
                project=project,
                dataset=dataset,
                extract_table_prefix=extract_table_prefix,
                bin_number=bin_num,
                output_path=output_path,
                extract_metadata=extract_metadata,
                obs_columns=obs_columns,
            )
            futures.append(future)

        # Wait for all extractions to complete
        concurrency.wait(futures, return_when=concurrency.ALL_COMPLETED)

        # Check for errors
        for future in futures:
            future.result()  # Will raise any exceptions that occurred
