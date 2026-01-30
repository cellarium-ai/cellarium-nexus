"""SOMA data operation pipelines for Kubeflow."""

import typing as t

from kfp import dsl

from cellarium.nexus.workflows.kubeflow import components


@dsl.pipeline(
    name="nexus-pipelines-run-soma-randomized-extracts",
    description="Run parallel SOMA randomized extract operations",
)
def run_soma_randomized_extracts_pipeline(extract_configs: t.List[str]) -> None:
    """
    Run multiple SOMA randomized extract operations in parallel.

    :param extract_configs: List of paths to configuration files for soma_randomized_extract_job

    :raise RuntimeError: If any component fails
    """
    with dsl.ParallelFor(items=extract_configs, name="soma-randomized-extract-workers", parallelism=256) as item:
        components.soma_randomized_extract_job(config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-run-soma-grouped-extracts",
    description="Run parallel SOMA grouped extract operations",
)
def run_soma_grouped_extracts_pipeline(extract_configs: t.List[str]) -> None:
    """
    Run multiple SOMA grouped extract operations in parallel.

    :param extract_configs: List of paths to configuration files for soma_grouped_extract_job

    :raise RuntimeError: If any component fails
    """
    with dsl.ParallelFor(items=extract_configs, name="soma-grouped-extract-workers", parallelism=256) as item:
        components.soma_grouped_extract_job(config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-soma-randomized-extract-data",
    description="Extract data from SOMA experiment using randomized extraction",
)
def soma_randomized_extract_data_pipeline(extract_configs: t.List[str], mark_finished_config: str) -> None:
    """
    Extract data from SOMA experiment using randomized extraction with shuffle.

    The extract plan and curriculum registration must already be created before
    this pipeline is invoked. This pipeline runs randomized extraction workers
    and marks the curriculum as finished.

    :param extract_configs: List of paths to configuration files for soma_randomized_extract_job
    :param mark_finished_config: Path to configuration file for mark_soma_curriculum_as_finished_job

    :raise RuntimeError: If any component fails
    """
    extract_op = run_soma_randomized_extracts_pipeline(extract_configs=extract_configs)

    mark_finished_op = components.mark_soma_curriculum_as_finished_job(config_path=mark_finished_config)
    mark_finished_op.after(extract_op)


@dsl.pipeline(
    name="nexus-pipelines-soma-grouped-extract-data",
    description="Extract data from SOMA experiment using grouped extraction",
)
def soma_grouped_extract_data_pipeline(extract_configs: t.List[str], mark_finished_config: str) -> None:
    """
    Extract data from SOMA experiment using grouped extraction.

    The extract plan and curriculum registration must already be created before
    this pipeline is invoked. This pipeline runs grouped extraction workers
    and marks the curriculum as finished.

    :param extract_configs: List of paths to configuration files for soma_grouped_extract_job
    :param mark_finished_config: Path to configuration file for mark_soma_curriculum_as_finished_job

    :raise RuntimeError: If any component fails
    """
    extract_op = run_soma_grouped_extracts_pipeline(extract_configs=extract_configs)

    mark_finished_op = components.mark_soma_curriculum_as_finished_job(config_path=mark_finished_config)
    mark_finished_op.after(extract_op)


@dsl.pipeline(
    name="nexus-pipelines-run-soma-validate-sanitize",
    description="Run parallel SOMA validate and sanitize operations",
)
def run_soma_validate_sanitize_pipeline(validate_sanitize_configs: t.List[str]) -> None:
    """
    Validate and sanitize multiple SOMA ingest inputs in parallel.

    :param validate_sanitize_configs: List of paths to SomaValidateSanitizeConfig files

    :raise RuntimeError: If any component fails
    """
    with dsl.ParallelFor(
        items=validate_sanitize_configs, name="soma-validate-sanitize-workers", parallelism=64
    ) as item:
        components.soma_validate_sanitize_job(config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-soma-validate-sanitize-and-ingest",
    description="Validate, sanitize, and ingest SOMA data",
)
def soma_validate_sanitize_and_ingest_pipeline(
    validate_sanitize_configs: t.List[str], ingest_plan_config: str, ingest_partition_configs: t.List[str]
) -> None:
    """
    Validate and sanitize inputs, prepare ingest plan, and ingest partitions.

    :param validate_sanitize_configs: List of paths to SomaValidateSanitizeConfig files
    :param ingest_plan_config: Path to SomaIngestPlanConfig file
    :param ingest_partition_configs: List of paths to SomaIngestPartitionConfig files

    :raise RuntimeError: If any component fails
    """
    validate_op = run_soma_validate_sanitize_pipeline(validate_sanitize_configs=validate_sanitize_configs)
    plan_op = components.soma_prepare_ingest_plan_job(config_path=ingest_plan_config)
    plan_op.after(validate_op)

    with dsl.ParallelFor(items=ingest_partition_configs, name="soma-ingest-partition-workers", parallelism=128) as item:
        ingest_op = components.soma_ingest_partition_job(config_path=item)
        ingest_op.after(plan_op)


@dsl.pipeline(
    name="nexus-pipelines-soma-ingest",
    description="Prepare ingest plan and ingest SOMA data",
)
def soma_ingest_pipeline(ingest_plan_config: str, ingest_partition_configs: t.List[str]) -> None:
    """
    Prepare ingest plan and ingest partitions.

    :param ingest_plan_config: Path to SomaIngestPlanConfig file
    :param ingest_partition_configs: List of paths to SomaIngestPartitionConfig files

    :raise RuntimeError: If any component fails
    """
    plan_op = components.soma_prepare_ingest_plan_job(config_path=ingest_plan_config)

    with dsl.ParallelFor(items=ingest_partition_configs, name="soma-ingest-partition-workers", parallelism=128) as item:
        ingest_op = components.soma_ingest_partition_job(config_path=item)
        ingest_op.after(plan_op)
