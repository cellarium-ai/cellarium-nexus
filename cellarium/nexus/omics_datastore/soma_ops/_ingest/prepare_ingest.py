import logging

import tiledbsoma.io
from anndata import AnnData
from tiledbsoma import DoesNotExistError

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import sanitize_for_ingest, validate_for_ingest
from cellarium.nexus.shared.schemas.omics_datastore import IngestSchema

logger = logging.getLogger(__name__)


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
    :param ingest_schema: Schema to validate the AnnData against.
    """
    logger.info("Validating AnnData for ingest")
    validate_for_ingest(adata=adata, schema=ingest_schema)
    logger.info("Sanitizing AnnData for ingest")
    sanitize_for_ingest(adata=adata)


def prepare_for_ingest(
    *,
    experiment_uri: str,
    h5ad_paths: list[str],
    measurement_name: str,
) -> tiledbsoma.io.ExperimentAmbientLabelMapping:
    """
    Prepare the SOMA experiment for ingestion.

    Create a schema-only experiment if it does not exist, register h5ad files,
    and prepare the experiment for append workflow.

    :param experiment_uri: URI of the SOMA experiment.
    :param h5ad_paths: List of paths to h5ad files to register.
    :param measurement_name: Name of the measurement.

    :returns: Experiment mapping for the registered h5ad files.
    """
    if not _soma_experiment_exists(soma_uri=experiment_uri):
        logger.info(f"Creating schema-only SOMA experiment at URI: {experiment_uri}")
        tiledbsoma.io.from_h5ad(
            experiment_uri=experiment_uri,
            input_path=h5ad_paths[0],
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


def ingest_h5ads(
    *,
    experiment_uri: str,
    h5ad_paths: list[str],
    measurement_name: str,
    registration_mapping: tiledbsoma.io.ExperimentAmbientLabelMapping,
) -> None:
    """
    Ingest h5ad files into a SOMA experiment.

    :param experiment_uri: URI of the SOMA experiment.
    :param h5ad_paths: List of paths to h5ad files to ingest.
    :param measurement_name: Name of the measurement.
    :param registration_mapping: Registration mapping from prepare_for_ingest.
    """
    logger.info(f"Starting ingestion of {len(h5ad_paths)} h5ad files")
    for idx, h5ad_path in enumerate(h5ad_paths, start=1):
        logger.info(f"Ingesting h5ad file {idx}/{len(h5ad_paths)}: {h5ad_path}")
        tiledbsoma.io.from_h5ad(
            experiment_uri=experiment_uri,
            input_path=h5ad_path,
            measurement_name=measurement_name,
            obs_id_name="obs_id",
            var_id_name="var_id",
            ingest_mode="write",
            registration_mapping=registration_mapping,
        )
    logger.info("Completed ingestion of all h5ad files")
