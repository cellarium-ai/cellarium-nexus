import json
import logging
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Iterator

import fastavro
from django.conf import settings
from django.db import transaction
from nexus.omics_datastore.bq_ops import constants

from cellarium.nexus.backend.cell_management import models as cell_models
from cellarium.nexus.backend.ingest_management import models as ingest_models
from cellarium.nexus.shared import utils

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def parse_metadata_extra(metadata_extra: str | None) -> dict:
    """
    Parse metadata_extra JSON string into a dictionary.

    :param metadata_extra: JSON string from Avro file

    :return: Parsed dictionary or empty dict if input is None or invalid
    """
    if not metadata_extra:
        return {}

    try:
        return json.loads(metadata_extra)
    except (json.JSONDecodeError, TypeError):
        logger.warning(f"Failed to parse metadata_extra: {metadata_extra}. Using empty dict instead.")
        return {}


def read_avro_records(file_path: str | Path) -> Iterator[dict]:
    """
    Read records from an Avro file.

    :param file_path: Path to the Avro file

    :raise FileNotFoundError: if file doesn't exist
    :raise Exception: for other errors

    :return: Iterator of records from the Avro file
    """
    try:
        with open(file_path, "rb") as avro_file:
            reader = fastavro.reader(avro_file)
            yield from reader
    except FileNotFoundError:
        logger.error(f"Avro file not found: {file_path}")
        raise
    except Exception as e:
        logger.error(f"Error reading Avro file {file_path}: {e}")
        raise


def fetch_cell_info(stage_dir: str, file_name: str) -> Sequence[cell_models.CellInfo]:
    """
    Fetch CellInfo records from an Avro file in GCS.

    :param stage_dir: Staging directory path containing the files
    :param file_name: Name of the Avro file

    :raise storage.exceptions.NotFound: if file doesn't exist in bucket
    :raise Exception: for other errors

    :return: A list of populated CellInfo models
    """
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            local_path = Path(temp_dir) / file_name

            source_blob_name = f"{stage_dir}/{file_name}"
            utils.gcp.download_file_from_bucket(
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                source_blob_name=source_blob_name,
                destination_file_name=local_path,
            )

            cell_info_list = []
            ingest_cache = {}

            for record in read_avro_records(local_path):
                ingest_id = record["ingest_id"]

                # Cache IngestInfo lookups to avoid repeated DB hits
                if ingest_id not in ingest_cache:
                    ingest_cache[ingest_id] = ingest_models.IngestInfo.objects.only("id", "bigquery_dataset_id").get(
                        id=ingest_id
                    )

                ingest = ingest_cache[ingest_id]

                cell_info = cell_models.CellInfo(
                    id=record["id"],
                    original_id=record["original_id"],
                    ingest_id=ingest.id,
                    bigquery_dataset_id=ingest.bigquery_dataset_id,
                    metadata_extra=parse_metadata_extra(record["metadata_extra"]),
                    donor_id=record.get("donor_id"),
                    cell_type=record["cell_type"],
                    assay=record.get("assay"),
                    development_stage=record.get("development_stage"),
                    tissue=record.get("tissue"),
                    disease=record.get("disease"),
                    organism=record.get("organism"),
                    self_reported_ethnicity=record.get("self_reported_ethnicity"),
                    sex=record.get("sex"),
                    suspension_type=record.get("suspension_type"),
                    total_mrna_umis=record.get("total_mrna_umis"),
                    cell_type_ontology_term_id=record.get("cell_type_ontology_term_id"),
                    assay_ontology_term_id=record.get("assay_ontology_term_id"),
                    development_stage_ontology_term_id=record.get("development_stage_ontology_term_id"),
                    tissue_ontology_term_id=record.get("tissue_ontology_term_id"),
                    disease_ontology_term_id=record.get("disease_ontology_term_id"),
                    organism_ontology_term_id=record.get("organism_ontology_term_id"),
                    self_reported_ethnicity_ontology_term_id=record.get("self_reported_ethnicity_ontology_term_id"),
                    sex_ontology_term_id=record.get("sex_ontology_term_id"),
                    tag=record.get("tag"),
                )
                cell_info_list.append(cell_info)

            logger.info(f"Fetched {len(cell_info_list)} CellInfo records from {source_blob_name}")
            return cell_info_list

    except Exception as e:
        logger.error(f"Error processing CellInfo records from {stage_dir}/{file_name}: {e}")
        raise


def fetch_feature_info(stage_dir: str, file_name: str) -> Sequence[cell_models.CellFeatureInfo]:
    """
    Fetch FeatureInfo records from an Avro file in GCS.

    :param stage_dir: Staging directory path containing the files
    :param file_name: Name of the Avro file

    :raise storage.exceptions.NotFound: if file doesn't exist in bucket
    :raise Exception: for other errors

    :return: A list of populated FeatureInfo models
    """
    try:
        # Create a temporary directory that will be automatically cleaned up
        with tempfile.TemporaryDirectory() as temp_dir:
            local_path = Path(temp_dir) / file_name

            # Download file from GCS using our utility
            source_blob_name = f"{stage_dir}/{file_name}"
            utils.gcp.download_file_from_bucket(
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                source_blob_name=source_blob_name,
                destination_file_name=local_path,
            )

            feature_info_list = []
            for record in read_avro_records(local_path):
                feature_info = cell_models.CellFeatureInfo(
                    id=record["id"],
                    ensemble_id=record["ensemble_id"],
                    symbol=record["symbol"],
                    biotype=record.get("biotype"),
                    is_filtered=record.get("is_filtered"),
                    reference=record["reference"],
                    ingest_id=record["ingest_id"],
                    metadata_extra=parse_metadata_extra(record["metadata_extra"]),
                    tag=record.get("tag"),
                )
                feature_info_list.append(feature_info)

            logger.info(f"Fetched {len(feature_info_list)} FeatureInfo records from {source_blob_name}")
            return feature_info_list

    except Exception as e:
        logger.error(f"Error processing FeatureInfo records from {stage_dir}/{file_name}: {e}")
        raise


def ingest_files(stage_dir: str) -> tuple[int, int]:
    """
    Ingest CellInfo and FeatureInfo from Avro files in a single transaction.

    :param stage_dir: Base staging directory path

    :raise storage.exceptions.NotFound: if files don't exist in bucket
    :raise Exception: for other errors

    :return: Tuple of (cell_info_count, feature_info_count)
    """
    with transaction.atomic():
        # Fetch records
        cell_infos = fetch_cell_info(
            stage_dir=stage_dir,
            file_name=constants.INGEST_CELL_INFO_FILE_NAME,
        )
        feature_infos = fetch_feature_info(
            stage_dir=stage_dir,
            file_name=constants.INGEST_FEATURE_INFO_FILE_NAME,
        )

        # Bulk insert into database
        cell_models.CellInfo.objects.bulk_create(objs=cell_infos, batch_size=settings.INGEST_BATCH_SIZE)
        cell_models.CellFeatureInfo.objects.bulk_create(objs=feature_infos, batch_size=settings.INGEST_BATCH_SIZE)

        cell_info_count = len(cell_infos)
        feature_info_count = len(feature_infos)

        logger.info(
            f"Successfully ingested {cell_info_count} Cell Info objects and {feature_info_count} Feature Info objects"
        )
        return cell_info_count, feature_info_count
