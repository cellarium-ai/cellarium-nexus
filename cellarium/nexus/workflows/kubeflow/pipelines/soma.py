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
def soma_extract_data_pipeline(prepare_extract_config: str, extract_configs: t.List[str]) -> None:
    """
    Extract data from SOMA experiment into AnnData files.

    :param prepare_extract_config: Path to configuration file for prepare_soma_extract_job
    :param extract_configs: List of paths to configuration files for soma_extract_job

    :raise RuntimeError: If any component fails
    """
    prepare_op = components.prepare_soma_extract_job(config_path=prepare_extract_config)

    extract_op = run_soma_extracts_pipeline(extract_configs=extract_configs)
    extract_op.after(prepare_op)

    mark_finished_op = components.mark_soma_curriculum_as_finished_job(config_path=prepare_extract_config)
    mark_finished_op.after(extract_op)
