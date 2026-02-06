from anndata import AnnData

from cellarium.nexus.omics_datastore.soma_ops._ingest.ingest import (
    ingest_h5ads_partition,
    prepare_ingest_plan,
    validate_and_sanitize_for_ingest,
)
from cellarium.nexus.shared.schemas.omics_datastore import IngestPlanMetadata, IngestSchema


class TileDBSOMAIngestor:

    def validate_and_sanitize_for_ingest(self, *, adata: AnnData, ingest_schema: IngestSchema) -> None:
        """
        Validate and sanitize an AnnData object for ingestion.

        :param adata: The AnnData object to validate and sanitize.
        :param ingest_schema: Schema to validate the AnnData against.
        """
        validate_and_sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)

    def prepare_ingest_plan(
        self,
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
        return prepare_ingest_plan(
            experiment_uri=experiment_uri,
            h5ad_paths=h5ad_paths,
            measurement_name=measurement_name,
            ingest_schema=ingest_schema,
            ingest_batch_size=ingest_batch_size,
            first_adata=first_adata,
        )

    def ingest_h5ads_partition(
        self,
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
        :param local_h5ad_paths: List of local file paths to h5ad files
            (pre-downloaded and pre-sliced by coordinator).

        :raises ValueError: If local_h5ad_paths is empty.
        """
        ingest_h5ads_partition(
            ingest_plan=ingest_plan,
            local_h5ad_paths=local_h5ad_paths,
        )
