"""SOMA data operation pipelines for Kubeflow."""

import typing as t

from kfp import dsl

from cellarium.nexus.workflows.kubeflow import components


@dsl.pipeline(
    name="nexus-pipelines-run-soma-extracts",
    description="Run parallel SOMA extract operations",
)
def run_soma_extracts_pipeline(extract_configs: t.List[str]) -> None:
    """
    Run multiple SOMA extract operations in parallel.

    :param extract_configs: List of paths to configuration files for soma_extract_job

    :raise RuntimeError: If any component fails
    """
    with dsl.ParallelFor(items=extract_configs, name="soma-extract-workers", parallelism=128) as item:
        components.soma_extract_job(config_path=item)


@dsl.pipeline(
    name="nexus-pipelines-soma-extract-data",
    description="Extract data from SOMA experiment",
)
def soma_extract_data_pipeline(extract_configs: t.List[str], mark_finished_config: str) -> None:
    """
    Extract data from SOMA experiment into AnnData files.

    The extract plan and curriculum registration must already be created before
    this pipeline is invoked (e.g., by the admin flow). This pipeline only
    runs extraction workers and marks the curriculum as finished.

    :param extract_configs: List of paths to configuration files for soma_extract_job

    :raise RuntimeError: If any component fails
    """
    extract_op = run_soma_extracts_pipeline(extract_configs=extract_configs)

    mark_finished_op = components.mark_soma_curriculum_as_finished_job(config_path=mark_finished_config)
    mark_finished_op.after(extract_op)
