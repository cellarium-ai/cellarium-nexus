"""
Utility functions for SOMA extract workflow submission.
"""

import math

from django.conf import settings

from cellarium.nexus.backend.cell_management import models
from cellarium.nexus.coordinator import SomaDataOpsCoordinator
from cellarium.nexus.omics_datastore.soma_ops import TileDBSOMADataOperator
from cellarium.nexus.shared.schemas import SomaCurriculumMetadata, component_configs
from cellarium.nexus.shared.utils import workflows_configs


def get_total_cells_in_soma(omics_dataset: models.OmicsDataset, filters: dict | None = None) -> int:
    """
    Get total number of cells in SOMA experiment that match the given filters.

    :param omics_dataset: Omics dataset with TileDB SOMA backend
    :param filters: Optional dictionary of filter statements to apply

    :raise ValueError: If omics_dataset has no URI configured
    :raise SomaReadError: If SOMA read operation fails

    :return: Total number of cells matching the filters
    """
    if not omics_dataset.uri:
        raise ValueError(f"TileDB SOMA dataset '{omics_dataset.name}' has no URI configured")

    soma_operator = TileDBSOMADataOperator(experiment_uri=omics_dataset.uri)
    return soma_operator.count_cells(filter_statements=filters)


def _get_extract_paths(name: str) -> tuple[str, str]:
    """
    Build extract bucket and extract_metadata paths for a SOMA extract.

    :param name: Name for the extract curriculum

    :raise ValueError: If name is empty

    :return: Tuple of extract_bucket_path and extract_metadata_path
    """
    if not name:
        raise ValueError("Extract name must be provided")

    extract_bucket_path = f"{settings.BACKEND_PIPELINE_DIR}/{settings.PIPELINE_DATA_EXTRACTS_DIR}/{name}"
    extract_metadata_path = f"{extract_bucket_path}/extract_metadata.json"

    return extract_bucket_path, extract_metadata_path


def compose_soma_extract_configs(
    name: str,
    omics_dataset: models.OmicsDataset,
    curriculum_metadata: SomaCurriculumMetadata,
    extract_metadata_path: str,
    extract_bucket_path: str,
    output_format: str = "h5ad",
    max_workers_extract: int | None = None,
    max_workers_shuffle: int | None = None,
) -> list[component_configs.SomaOpsExtract]:
    """
    Compose SOMA extract configs for the Kubeflow pipeline.

    Create one config per worker, where number of workers is determined by
    ranges per worker setting. Each worker gets a consecutive partition index.

    :param name: Name for the extract curriculum
    :param omics_dataset: Omics dataset with TileDB SOMA backend
    :param curriculum_metadata: SOMA extract metadata
    :param extract_metadata_path: Bucket path where worker can obtain the extract metadata
    :param extract_bucket_path: Bucket path where to save output chunks
    :param output_format: Output format - "h5ad" or "zarr"
    :param max_workers_extract: Maximum parallel workers for extraction
    :param max_workers_shuffle: Maximum parallel workers for shuffling

    :raise ValueError: If omics_dataset has no URI

    :return: List of extract configs
    """
    if not omics_dataset.uri:
        raise ValueError(f"TileDB SOMA dataset '{omics_dataset.name}' has no URI configured")

    num_ranges = len(curriculum_metadata.id_ranges)
    max_ranges_per_worker = settings.TILEDB_SOMA_RANGES_PER_WORKER
    num_workers = math.ceil(num_ranges / max_ranges_per_worker)

    return [
        component_configs.SomaOpsExtract(
            extract_name=name,
            experiment_uri=omics_dataset.uri,
            nexus_backend_api_url=settings.SITE_URL,
            bucket_name=settings.BUCKET_NAME_PRIVATE,
            extract_metadata_path=extract_metadata_path,
            extract_bucket_path=extract_bucket_path,
            curriculum_partition_index=i,
            curriculum_partitions_num=num_workers,
            output_format=output_format,
            max_workers_extract=max_workers_extract,
            max_workers_shuffle=max_workers_shuffle,
        )
        for i in range(num_workers)
    ]


def compose_and_dump_soma_configs(
    name: str,
    creator_id: int,
    omics_dataset: models.OmicsDataset,
    range_size: int,
    output_chunk_size: int,
    feature_schema: models.FeatureSchema | None = None,
    filters: dict | None = None,
    shuffle_ranges: bool = True,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    x_layer: str = "raw",
    output_format: str = "h5ad",
    max_workers_extract: int | None = None,
    max_workers_shuffle: int | None = None,
    prepare: bool = True,
) -> list[str]:
    """
    Compose SOMA extract configs, create the extract plan, and dump configs to GCS bucket.

    :param name: Name for the extract curriculum
    :param creator_id: ID of the user who initiated the extract
    :param omics_dataset: Omics dataset with TileDB SOMA backend
    :param range_size: Target number of cells per range
    :param output_chunk_size: Target cells per output chunk (for shuffling)
    :param feature_schema: Feature schema used to derive feature IDs for filtering
    :param filters: Optional dictionary of filter statements to apply
    :param shuffle_ranges: Whether to shuffle the joinid ranges
    :param obs_columns: Optional obs columns to include in output files
    :param var_columns: Optional var columns to include in output files
    :param x_layer: Name of the SOMA X layer to read counts from
    :param output_format: Output format - "h5ad" or "zarr"
    :param max_workers_extract: Maximum parallel workers for extraction
    :param max_workers_shuffle: Maximum parallel workers for shuffling

    :raise ValueError: If range_size <= 0 or output_chunk_size <= 0 or omics_dataset has no URI
    :raise exceptions.ZeroCellsReturnedError: If no cells match the filters
    :raise IOError: If there's an error writing configs to GCS
    :raise exceptions.SomaExtractError: If planning fails
    :raise HTTPError: If backend API calls fail

    :return: List of extract config paths
    """
    var_filter_values = [feature.ensemble_id for feature in feature_schema.features.all()] if feature_schema else None
    var_filter_column = "feature_id" if var_filter_values is not None else None

    extract_bucket_path, extract_metadata_path = _get_extract_paths(name=name)
    coordinator = SomaDataOpsCoordinator(
        experiment_uri=omics_dataset.uri,
        nexus_backend_api_url=settings.SITE_URL,
        bucket_name=settings.BUCKET_NAME_PRIVATE,
    )
    curriculum_metadata = coordinator.prepare_soma_extract(
        extract_name=name,
        creator_id=creator_id,
        curriculum_metadata_path=extract_metadata_path,
        range_size=range_size,
        output_chunk_size=output_chunk_size,
        filters=filters,
        shuffle_ranges=shuffle_ranges,
        var_filter_column=var_filter_column,
        var_filter_values=var_filter_values,
        obs_columns=obs_columns,
        var_columns=var_columns,
        x_layer=x_layer,
    )

    extract_configs = compose_soma_extract_configs(
        name=name,
        omics_dataset=omics_dataset,
        extract_bucket_path=extract_bucket_path,
        curriculum_metadata=curriculum_metadata,
        extract_metadata_path=extract_metadata_path,
        output_format=output_format,
        max_workers_extract=max_workers_extract,
        max_workers_shuffle=max_workers_shuffle,
    )

    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/{settings.PIPELINE_CONFIGS_DIR}"

    extract_config_paths = workflows_configs.dump_configs_to_bucket(
        configs=extract_configs, bucket_path=configs_stage_dir
    )

    return extract_config_paths


def submit_soma_extract_pipeline(
    name: str,
    creator_id: int,
    omics_dataset: models.OmicsDataset,
    range_size: int,
    output_chunk_size: int,
    feature_schema: models.FeatureSchema | None = None,
    filters: dict | None = None,
    shuffle_ranges: bool = True,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    x_layer: str = "raw",
    output_format: str = "h5ad",
    max_workers_extract: int | None = None,
    max_workers_shuffle: int | None = None,
) -> str:
    """
    Submit a Kubeflow pipeline for SOMA data extraction.

    Generate SOMA extract configs, save them to GCS, and submit the pipeline for execution.

    :param name: Name for the extract curriculum
    :param creator_id: ID of the user who initiated the extract
    :param omics_dataset: Omics dataset with TileDB SOMA backend
    :param range_size: Target number of cells per range
    :param output_chunk_size: Target cells per output chunk (for shuffling)
    :param feature_schema: Feature schema used to derive feature IDs for filtering
    :param var_filter_column: Column to use for filtering the features
    :param var_filter_values: Values to match in the feature filter column
    :param filters: Optional dictionary of filter statements to apply
    :param shuffle_ranges: Whether to shuffle the joinid ranges
    :param obs_columns: Optional obs columns to include in output files
    :param var_columns: Optional var columns to include in output files
    :param x_layer: Name of the SOMA X layer to read counts from
    :param output_format: Output format - "h5ad" or "zarr"
    :param max_workers_extract: Maximum parallel workers for extraction
    :param max_workers_shuffle: Maximum parallel workers for shuffling

    :raise ValueError: If range_size <= 0 or output_chunk_size <= 0 or omics_dataset has no URI
    :raise exceptions.ZeroCellsReturnedError: If no cells match the filters
    :raise IOError: If there's an error writing configs to GCS
    :raise google.api_core.exceptions.GoogleAPIError: If pipeline submission fails

    :return: URL to the Vertex AI Pipeline dashboard for the submitted job
    """
    from cellarium.nexus.workflows.kubeflow.pipelines import soma_extract_data_pipeline
    from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline

    extract_config_paths = compose_and_dump_soma_configs(
        name=name,
        creator_id=creator_id,
        omics_dataset=omics_dataset,
        range_size=range_size,
        output_chunk_size=output_chunk_size,
        feature_schema=feature_schema,
        filters=filters,
        shuffle_ranges=shuffle_ranges,
        obs_columns=obs_columns,
        var_columns=var_columns,
        x_layer=x_layer,
        output_format=output_format,
        max_workers_extract=max_workers_extract,
        max_workers_shuffle=max_workers_shuffle,
    )

    return submit_pipeline(
        pipeline_component=soma_extract_data_pipeline,
        display_name=f"Nexus SOMA Extract - {name}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={
            "extract_configs": extract_config_paths,
            "mark_finished_config": extract_config_paths[0],
        },
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
        labels={"application": settings.GCP_APPLICATION_BILLING_LABEL, "method": "soma_extract"},
    )


def run_soma_extract_pipeline_locally(
    name: str,
    creator_id: int,
    omics_dataset: models.OmicsDataset,
    range_size: int,
    output_chunk_size: int,
    var_filter_column: str,
    feature_schema: models.FeatureSchema,
    filters: dict | None = None,
    shuffle_ranges: bool = True,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    x_layer: str = "raw",
    output_format: str = "h5ad",
    max_workers_extract: int | None = None,
    max_workers_shuffle: int | None = None,
) -> None:
    """
    Run a SOMA extract pipeline locally for testing and development.

    Execute the SOMA extract workflow by directly calling the coordinator methods
    instead of submitting to Vertex AI.

    :param name: Name for the extract curriculum
    :param creator_id: ID of the user who initiated the extract
    :param omics_dataset: Omics dataset with TileDB SOMA backend
    :param range_size: Target number of cells per range
    :param var_filter_column: Column to use for filtering the features (used with `var_filter_values`)
    :param feature_schema: Object of feature schema
    :param output_chunk_size: Target cells per output chunk (for shuffling)
    :param filters: Optional dictionary of filter statements to apply
    :param shuffle_ranges: Whether to shuffle the joinid ranges
    :param obs_columns: Optional obs columns to include in output files
    :param var_columns: Optional var columns to include in output files
    :param x_layer: Name of the SOMA X layer to read counts from
    :param output_format: Output format - "h5ad" or "zarr"
    :param max_workers_extract: Maximum parallel workers for extraction
    :param max_workers_shuffle: Maximum parallel workers for shuffling

    :raise ValueError: If range_size <= 0 or output_chunk_size <= 0 or omics_dataset has no URI
    :raise exceptions.ZeroCellsReturnedError: If no cells match the filters
    :raise SomaExtractError: If extraction fails
    """

    if not omics_dataset.uri:
        raise ValueError(f"TileDB SOMA dataset '{omics_dataset.name}' has no URI configured")

    extract_bucket_path = f"{settings.BACKEND_PIPELINE_DIR}/{settings.PIPELINE_DATA_EXTRACTS_DIR}/{name}"
    extract_metadata_path = f"{extract_bucket_path}/extract_metadata.json"

    coordinator = SomaDataOpsCoordinator(
        experiment_uri=omics_dataset.uri,
        nexus_backend_api_url=settings.SITE_URL,
        bucket_name=settings.BUCKET_NAME_PRIVATE,
    )

    # Step 1: Prepare extract (compute metadata, register curriculum)
    feature_ids = [feature.ensemble_id for idx, feature in enumerate(feature_schema.features.all())]
    coordinator.prepare_soma_extract(
        extract_name=name,
        creator_id=creator_id,
        curriculum_metadata_path=extract_metadata_path,
        range_size=range_size,
        output_chunk_size=output_chunk_size,
        filters=filters,
        shuffle_ranges=shuffle_ranges,
        var_filter_column=var_filter_column,
        var_filter_values=feature_ids,
        obs_columns=obs_columns,
        var_columns=var_columns,
        x_layer=x_layer,
    )

    # Step 2: Run extraction
    coordinator.run_soma_extract(
        extract_name=name,
        curriculum_metadata_path=extract_metadata_path,
        extract_bucket_path=extract_bucket_path,
        output_format=output_format,
        max_workers_extract=max_workers_extract,
        max_workers_shuffle=max_workers_shuffle,
    )

    # Step 3: Mark curriculum as finished
    coordinator.mark_soma_curriculum_as_finished(
        extract_name=name,
        curriculum_metadata_path=extract_metadata_path,
        extract_bucket_path=extract_bucket_path,
    )
