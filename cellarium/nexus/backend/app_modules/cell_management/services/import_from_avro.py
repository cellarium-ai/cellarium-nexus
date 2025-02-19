import json
import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Iterator

import fastavro
from django.conf import settings
from django.db import transaction
from nexus.backend.app_modules.cell_management import models
from nexus.omics_datastore import gc_utils
from nexus.omics_datastore.bq_ops import constants

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


def fetch_cell_info(stage_dir: str, ingest_uuid: str, file_name: str) -> Sequence[models.CellInfo]:
    """
    Fetch CellInfo records from an Avro file in GCS.

    :param stage_dir: Base staging directory path
    :param ingest_uuid: UUID of the ingest process
    :param file_name: Name of the Avro file

    :raise storage.exceptions.NotFound: if file doesn't exist in bucket
    :raise Exception: for other errors

    :return: A list of populated CellInfo models
    """
    try:
        # Create temporary directory
        temp_dir = Path("/tmp") / str(ingest_uuid)
        temp_dir.mkdir(parents=True, exist_ok=True)
        local_path = temp_dir / file_name

        # Download file from GCS using our utility
        source_blob_name = f"{stage_dir}/{ingest_uuid}/{file_name}"
        gc_utils.download_file_from_bucket(
            bucket_name=settings.BUCKET_NAME_PRIVATE,
            source_blob_name=source_blob_name,
            destination_file_name=local_path,
        )

        try:
            cell_info_list = []
            for record in read_avro_records(local_path):
                cell_info = models.CellInfo(
                    id=record["id"],
                    original_id=record["original_id"],
                    ingest_id=record["ingest_id"],
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

        finally:
            # Clean up temporary file
            local_path.unlink(missing_ok=True)
            temp_dir.rmdir()

    except Exception as e:
        logger.error(f"Error processing CellInfo records from {stage_dir}/{ingest_uuid}/{file_name}: {e}")
        raise


def fetch_feature_info(stage_dir: str, ingest_uuid: str, file_name: str) -> Sequence[models.CellFeatureInfo]:
    """
    Fetch FeatureInfo records from an Avro file in GCS.

    :param stage_dir: Base staging directory path
    :param ingest_uuid: UUID of the ingest process
    :param file_name: Name of the Avro file

    :raise storage.exceptions.NotFound: if file doesn't exist in bucket
    :raise Exception: for other errors

    :return: A list of populated FeatureInfo models
    """
    try:
        # Create temporary directory
        temp_dir = Path("/tmp") / str(ingest_uuid)
        temp_dir.mkdir(parents=True, exist_ok=True)
        local_path = temp_dir / file_name

        # Download file from GCS using our utility
        source_blob_name = f"{stage_dir}/{ingest_uuid}/{file_name}"
        gc_utils.download_file_from_bucket(
            bucket_name=settings.BUCKET_NAME_PRIVATE,
            source_blob_name=source_blob_name,
            destination_file_name=local_path,
        )

        try:
            feature_info_list = []
            for record in read_avro_records(local_path):
                feature_info = models.CellFeatureInfo(
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

        finally:
            # Clean up temporary file
            local_path.unlink(missing_ok=True)
            temp_dir.rmdir()

    except Exception as e:
        logger.error(f"Error processing FeatureInfo records from {stage_dir}/{ingest_uuid}/{file_name}: {e}")
        raise


def ingest_files(stage_dir: str, ingest: models.IngestInfo) -> tuple[int, int]:
    """
    Ingest CellInfo and FeatureInfo from Avro files in a single transaction.

    :param stage_dir: Base staging directory path
    :param ingest: IngestFileInfo instance

    :raise storage.exceptions.NotFound: if files don't exist in bucket
    :raise Exception: for other errors

    :return: Tuple of (cell_info_count, feature_info_count)
    """
    with transaction.atomic():
        # Fetch records
        cell_infos = fetch_cell_info(
            stage_dir=stage_dir, ingest_uuid=str(ingest.nexus_uuid), file_name=constants.INGEST_CELL_INFO_FILE_NAME
        )
        feature_infos = fetch_feature_info(
            stage_dir=stage_dir, ingest_uuid=str(ingest.nexus_uuid), file_name=constants.INGEST_FEATURE_INFO_FILE_NAME
        )

        # Bulk insert into database
        models.CellInfo.objects.bulk_create(cell_infos)
        models.CellFeatureInfo.objects.bulk_create(feature_infos)

        cell_info_count = len(cell_infos)
        feature_info_count = len(feature_infos)

        logger.info(
            f"Successfully ingested {cell_info_count} Cell Info objects and {feature_info_count} Feature Info objects"
        )
        return cell_info_count, feature_info_count
