"""Validation pipelines for Kubeflow."""

import typing as t

from kfp import dsl

from cellarium.nexus.workflows.kubeflow import components


@dsl.pipeline(
    name="nexus-pipelines-validate-ingest-files",
    description="Validate ingest AnnData files and report results",
)
def validate_anndata_pipeline(validation_configs: t.List[str]) -> None:
    """
    Validate multiple AnnData files and report validation results.

    :param validation_configs: GCS path to configuration file for validate_anndata_files_job

    :raise RuntimeError: If the validation component fails
    """
    with dsl.ParallelFor(items=validation_configs, name="validate-ingest-files-workers", parallelism=64) as item:
        components.validate_anndata_files_job(gcs_config_path=item)
