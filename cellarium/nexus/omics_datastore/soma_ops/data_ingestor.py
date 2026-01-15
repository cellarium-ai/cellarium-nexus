from anndata import AnnData

from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import sanitize_for_ingest, validate_for_ingest
from cellarium.nexus.shared.schemas.omics_datastore import IngestSchema


class TileDBSOMAIngestor:

    def validate_and_sanitize_for_ingest(self, *, adata: AnnData, ingest_schema: IngestSchema) -> None:
        """
        Validate and sanitize an AnnData object for ingestion.

        :param adata: The AnnData object to validate and sanitize.
        """
        validate_for_ingest(adata=adata, schema=ingest_schema)
        sanitize_for_ingest(adata=adata)
