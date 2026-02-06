import base64
import logging
import math
import pickle

import anndata
import tiledbsoma.io
from anndata import AnnData
from tiledbsoma import DoesNotExistError

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import (
    sanitize_first_adata_for_schema,
    sanitize_for_ingest,
    validate_for_ingest,
)
from cellarium.nexus.shared.schemas.omics_datastore import IngestPlanMetadata, IngestSchema

logger = logging.getLogger(__name__)


def serialize_registration_mapping(
    *,
    mapping: tiledbsoma.io.ExperimentAmbientLabelMapping,
) -> str:
    """
    Serialize an ExperimentAmbientLabelMapping to a base64-encoded string.

    :param mapping: The registration mapping to serialize.

    :returns: Base64-encoded pickle string.
    """
    return base64.b64encode(pickle.dumps(mapping)).decode("ascii")


def deserialize_registration_mapping(
    *,
    data: str,
) -> tiledbsoma.io.ExperimentAmbientLabelMapping:
    """
    Deserialize an ExperimentAmbientLabelMapping from a base64-encoded string.

    :param data: Base64-encoded pickle string.

    :returns: Reconstructed ExperimentAmbientLabelMapping.
    """
    return pickle.loads(base64.b64decode(data))


def _soma_experiment_exists(*, soma_uri: str) -> bool:
    """
    Check if a SOMA experiment exists at the given URI.

    :param soma_uri: URI of the SOMA experiment to check.

    :returns: True if the experiment exists, False otherwise.
    """
    try:
        tiledbsoma.Experiment.open(uri=soma_uri, mode="r").close()
        return True
    except DoesNotExistError:
        logger.debug(f"SOMA experiment does not exist at URI: {soma_uri}")
        return False


def validate_and_sanitize_for_ingest(*, adata: AnnData, ingest_schema: IngestSchema) -> None:
    """
    Validate and sanitize an AnnData object for ingestion.

    :param adata: AnnData object to validate and sanitize.
    :param ingest_schema: Schema to validate and sanitize the AnnData against.
    """
    logger.info("Validating AnnData for ingest")
    validate_for_ingest(adata=adata, schema=ingest_schema)
    logger.info("Sanitizing AnnData for ingest")
    sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)


def prepare_for_ingest(
    *,
    experiment_uri: str,
    h5ad_paths: list[str],
    measurement_name: str,
    first_adata: AnnData | None = None,
) -> tiledbsoma.io.ExperimentAmbientLabelMapping:
    """
    Prepare the SOMA experiment for ingestion.

    Create a schema-only experiment if it does not exist, register h5ad files,
    and prepare the experiment for append workflow.

    :param experiment_uri: URI of the SOMA experiment.
    :param h5ad_paths: List of paths to h5ad files to register.
    :param measurement_name: Name of the measurement.
    :param first_adata: AnnData object to use for schema creation. Required if
        the experiment does not exist. Should be pre-sanitized with the full
        feature schema using sanitize_first_adata_for_schema.

    :raises ValueError: If the experiment does not exist and first_adata is None.

    :returns: Experiment mapping for the registered h5ad files.
    """
    if not _soma_experiment_exists(soma_uri=experiment_uri):
        if first_adata is None:
            raise ValueError(
                "first_adata is required when the SOMA experiment does not exist. "
                "Provide an AnnData object sanitized with the full feature schema."
            )
        logger.info(f"Creating schema-only SOMA experiment at URI: {experiment_uri}")
        tiledbsoma.io.from_anndata(
            experiment_uri=experiment_uri,
            anndata=first_adata,
            measurement_name=measurement_name,
            obs_id_name="obs_id",
            var_id_name="var_id",
            ingest_mode="schema_only",
        )
    else:
        logger.info(f"SOMA experiment already exists at URI: {experiment_uri}")

    logger.info(f"Registering {len(h5ad_paths)} h5ad files for experiment")
    experiment_mapping = tiledbsoma.io.register_h5ads(
        experiment_uri=experiment_uri,
        h5ad_file_names=h5ad_paths,
        measurement_name=measurement_name,
        obs_field_name="obs_id",
        var_field_name="var_id",
        use_multiprocessing=True,
    )

    logger.info("Preparing experiment for append workflow")
    experiment_mapping.prepare_experiment(experiment_uri)
    return experiment_mapping


def prepare_ingest_plan(
    *,
    experiment_uri: str,
    h5ad_paths: list[str],
    measurement_name: str,
    ingest_schema: IngestSchema,
    ingest_batch_size: int,
    first_adata: AnnData | None = None,
) -> IngestPlanMetadata:
    """
    Prepare an ingest plan for partitioned SOMA ingestion.

    Create the SOMA experiment if needed, register all h5ad files, and compute
    partition information for parallel ingestion.

    :param experiment_uri: URI of the SOMA experiment.
    :param h5ad_paths: List of paths to h5ad files (GCS URIs or local paths).
    :param measurement_name: Name of the measurement.
    :param ingest_schema: Schema for validating AnnData objects during ingest.
    :param ingest_batch_size: Number of h5ad files to ingest per partition.
    :param first_adata: AnnData object to use for schema creation. Required if
        the experiment does not exist. Will be validated and sanitized with the
        full feature schema before use.

    :raises ValueError: If ingest_batch_size is not positive.
    :raises ValueError: If the experiment does not exist and first_adata is None.

    :returns: IngestPlanMetadata containing all info needed for partitioned ingest.
    """
    if ingest_batch_size <= 0:
        raise ValueError(f"ingest_batch_size must be positive, got {ingest_batch_size}")

    total_files = len(h5ad_paths)
    num_partitions = math.ceil(total_files / ingest_batch_size) if total_files > 0 else 0

    if num_partitions > 0:
        last_partition_size = total_files - (num_partitions - 1) * ingest_batch_size
    else:
        last_partition_size = 0

    logger.info(f"Preparing ingest plan for {total_files} files across {num_partitions} partitions")

    # Validate and sanitize first_adata if provided
    sanitized_first_adata = None
    if first_adata is not None:
        logger.info("Validating first AnnData for schema creation")
        validate_for_ingest(adata=first_adata, schema=ingest_schema)
        logger.info("Sanitizing first AnnData with full feature schema")
        sanitize_first_adata_for_schema(adata=first_adata, ingest_schema=ingest_schema)
        sanitized_first_adata = first_adata

    registration_mapping = prepare_for_ingest(
        experiment_uri=experiment_uri,
        h5ad_paths=h5ad_paths,
        measurement_name=measurement_name,
        first_adata=sanitized_first_adata,
    )

    serialized_mapping = serialize_registration_mapping(mapping=registration_mapping)

    return IngestPlanMetadata(
        experiment_uri=experiment_uri,
        source_h5ad_uris=h5ad_paths,
        measurement_name=measurement_name,
        total_files=total_files,
        ingest_batch_size=ingest_batch_size,
        num_partitions=num_partitions,
        last_partition_size=last_partition_size,
        ingest_schema=ingest_schema,
        registration_mapping_pickle=serialized_mapping,
    )


def ingest_h5ads_partition(
    *,
    ingest_plan: IngestPlanMetadata,
    local_h5ad_paths: list[str],
) -> None:
    """
    Ingest h5ad files into a SOMA experiment.

    Load each h5ad file, validate and sanitize it using the ingest schema,
    then ingest into SOMA.

    :param ingest_plan: The ingest plan metadata containing experiment info,
        schema, and serialized registration mapping.
    :param local_h5ad_paths: List of local file paths to h5ad files.
        These should be pre-downloaded by the coordinator from the
        corresponding URIs in `ingest_plan.source_h5ad_uris`.

    :raises ValueError: If local_h5ad_paths is empty.
    """
    if not local_h5ad_paths:
        logger.info("No files to process")
        return

    file_count = len(local_h5ad_paths)
    logger.info(f"Starting ingest: processing {file_count} files")

    registration_mapping = deserialize_registration_mapping(data=ingest_plan.registration_mapping_pickle)

    for idx, h5ad_path in enumerate(local_h5ad_paths, start=1):
        logger.info(f"Processing h5ad file {idx}/{file_count}: {h5ad_path}")

        adata = anndata.read_h5ad(filename=h5ad_path)
        validate_and_sanitize_for_ingest(adata=adata, ingest_schema=ingest_plan.ingest_schema)

        logger.info(f"Ingesting h5ad file {idx}/{file_count}: {h5ad_path}")
        tiledbsoma.io.from_anndata(
            experiment_uri=ingest_plan.experiment_uri,
            anndata=adata,
            measurement_name=ingest_plan.measurement_name,
            obs_id_name="obs_id",
            var_id_name="var_id",
            ingest_mode="write",
            registration_mapping=registration_mapping,
        )

    logger.info(f"Completed ingest: {file_count} files processed")
