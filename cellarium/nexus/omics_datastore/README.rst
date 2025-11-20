Nexus Omics Datastore
=====================
This module is responsible for data operations related to raw count matrices and other omics data in single-cell
analysis. The Omics Datastore provides multiple backend implementations:

- **bq_ops/**: BigQuery-based operations for data ingestion, extraction, and querying
- **soma_ops/**: TileDB SOMA-based operations for efficient data extraction (future)

The module includes scripts for preprocessing and postprocessing the data to harmonize and unify it.