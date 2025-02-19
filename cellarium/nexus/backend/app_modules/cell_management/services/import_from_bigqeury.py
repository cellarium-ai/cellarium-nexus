import json
import logging
from collections.abc import Sequence

from google.api_core.exceptions import BadRequest, NotFound
from google.cloud import bigquery

from cellarium.nexus.backend.app_modules.cell_management import models
from cellarium.nexus.scripts.bq_ops import constants

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def get_bigquery_client() -> bigquery.Client:
    """
    Initialize and return a BigQuery client.

    :raise: Exception if client initialization fails

    :return: An authenticated BigQuery client
    """
    try:
        client = bigquery.Client()
        logger.info("BigQuery client initialized successfully.")
        return client
    except Exception as e:
        logger.error(f"Failed to initialize BigQuery client: {e}")
        raise


def parse_metadata_extra(metadata_extra: str | None) -> dict:
    """Parse metadata_extra JSON string into a dictionary.

    :param metadata_extra: JSON string from BigQuery

    :return: Parsed dictionary or empty dict if input is None or invalid
    """
    if not metadata_extra:
        return {}

    try:
        return json.loads(metadata_extra)
    except (json.JSONDecodeError, TypeError):
        logger.warning(f"Failed to parse metadata_extra: {metadata_extra}. Using empty dict instead.")
        return {}


def fetch_cell_info(ingest_id: str, project: str, dataset: str) -> Sequence[models.CellInfo]:
    """
    Fetch CellInfo records from BigQuery filtered by the specified ingest_id.

    :param ingest_id: The ingest ID to filter the records
    :param project: GCP project ID
    :param dataset: BigQuery dataset name

    :raise NotFound: if table doesn't exist
    :raise BadRequest: if query is invalid
    :raise Exception: for other errors

    :return: A list of populated CellInfo models
    """
    client = get_bigquery_client()

    query = f"""
        SELECT *
        FROM `{project}.{dataset}.{constants.BQ_CELL_INFO_TABLE_NAME}`
        WHERE ingest_id = @ingest_id
    """

    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("ingest_id", "STRING", ingest_id)]
    )

    try:
        query_job = client.query(query, job_config=job_config)
        results = query_job.result()
        cell_info_list = []

        for row in results:
            cell_info = models.CellInfo(
                id=row.id,
                original_id=row.original_id,
                ingest_id=row.ingest_id,
                metadata_extra=parse_metadata_extra(row.metadata_extra),
                donor_id=row.donor_id,
                cell_type=row.cell_type,
                assay=row.assay,
                development_stage=row.development_stage,
                tissue=row.tissue,
                disease=row.disease,
                organism=row.organism,
                self_reported_ethnicity=row.self_reported_ethnicity,
                sex=row.sex,
                suspension_type=row.suspension_type,
                total_mrna_umis=row.total_mrna_umis,
                cell_type_ontology_term_id=row.cell_type_ontology_term_id,
                assay_ontology_term_id=row.assay_ontology_term_id,
                development_stage_ontology_term_id=row.development_stage_ontology_term_id,
                tissue_ontology_term_id=row.tissue_ontology_term_id,
                disease_ontology_term_id=row.disease_ontology_term_id,
                organism_ontology_term_id=row.organism_ontology_term_id,
                self_reported_ethnicity_ontology_term_id=row.self_reported_ethnicity_ontology_term_id,
                sex_ontology_term_id=row.sex_ontology_term_id,
            )
            cell_info_list.append(cell_info)

        logger.info(f"Fetched {len(cell_info_list)} CellInfo records for ingest_id: {ingest_id}")
        return cell_info_list

    except NotFound:
        logger.error("The specified table does not exist.")
        raise
    except BadRequest as e:
        logger.error(f"Bad request error: {e}")
        raise
    except Exception as e:
        logger.error(f"An error occurred while fetching CellInfo: {e}")
        raise


def fetch_feature_info(ingest_id: str, project: str, dataset: str) -> Sequence[models.CellFeatureInfo]:
    """
    Fetch FeatureInfo records from BigQuery filtered by the specified ingest_id.

    :param ingest_id: The ingest ID to filter the records
    :param project: GCP project ID
    :param dataset: BigQuery dataset name

    :raise NotFound: if table doesn't exist
    :raise BadRequest: if query is invalid
    :raise Exception: for other errors

    :return: A list of populated FeatureInfo models
    """
    client = get_bigquery_client()

    query = f"""
        SELECT *
        FROM `{project}.{dataset}.{constants.BQ_FEATURE_INFO_TABLE_NAME}`
        WHERE ingest_id = @ingest_id
    """

    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("ingest_id", "STRING", ingest_id)]
    )

    try:
        query_job = client.query(query, job_config=job_config)
        results = query_job.result()
        feature_info_list = []

        for row in results:
            feature_info = models.CellFeatureInfo(
                id=row.id,
                ensemble_id=row.ensemble_id,
                symbol=row.symbol,
                biotype=row.biotype,
                is_filtered=row.is_filtered,
                reference=row.reference,
                ingest_id=row.ingest_id,
                metadata_extra=parse_metadata_extra(row.metadata_extra),
            )
            feature_info_list.append(feature_info)

        logger.info(f"Fetched {len(feature_info_list)} FeatureInfo records for ingest_id: {ingest_id}")
        return feature_info_list

    except NotFound:
        logger.error("The specified table does not exist.")
        raise
    except BadRequest as e:
        logger.error(f"Bad request error: {e}")
        raise
    except Exception as e:
        logger.error(f"An error occurred while fetching FeatureInfo: {e}")
        raise
