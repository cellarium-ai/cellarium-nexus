import datetime
import os
import secrets
from typing import Any

import pandas as pd
from django.conf import settings

from cellarium.nexus.backend.cell_management.models import OmicsDataset
from cellarium.nexus.shared import schemas, utils


def submit_ingest_pipeline(
    df_ingest_file_info: pd.DataFrame,
    omics_dataset: OmicsDataset,
    column_mapping: dict[str, Any],
    validation_methods: list[str] | None = None,
) -> str:
    """
    Submit a Kubeflow pipeline for data ingestion with the provided configurations.

    Create task configs for each file in the dataframe, save them to GCS, and submit the pipeline.

    :param df_ingest_file_info: DataFrame containing information about files to ingest, must include gcs_file_path column
    :param omics_dataset: Omics dataset where data will be ingested
    :param column_mapping: Dictionary mapping input columns to schema columns
    :param validation_methods: List of validation method names to apply to each file, if provided.

    :raise IOError: If there's an error writing configs to GCS
    :raise google.api_core.exceptions.GoogleAPIError: If pipeline submission fails

    :return: URL to the Vertex AI Pipeline dashboard for the submitted job
    """
    from cellarium.nexus.workflows.kubeflow.pipelines import ingest_data_to_bigquery_pipeline  # noqa
    from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline  # noqa

    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    base_stage_dir = f"{settings.BACKEND_PIPELINE_DIR}/data-ingests/{timestamp}_{secrets.token_hex(12)}"

    # Create list for combined task configs
    create_ingest_files_configs = []

    for i, row in df_ingest_file_info.iterrows():
        # Extract file name without extension from gcs_file_path
        file_path = row["gcs_file_path"]
        file_name = os.path.basename(file_path)
        file_name_without_ext = os.path.splitext(file_name)[0]

        # Create unique stage dir for this file
        stage_dir = f"{base_stage_dir}/{secrets.token_hex(8)}_{file_name_without_ext[:10]}"

        tag = row["tag"] if "tag" in df_ingest_file_info.columns else None

        # Create combined config for this task
        create_ingest_files_configs.append(
            schemas.component_configs.CreateIngestFilesConfig(
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=omics_dataset.name,
                input_file_path=file_path,
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                bucket_stage_dir=stage_dir,
                tag=tag,
                column_mapping=column_mapping,
                max_input_data_size=settings.INGEST_INPUT_FILE_MAX_SIZE,
                validation_methods=validation_methods,
                uns_keys_to_keep=settings.INGEST_UNS_KEYS_TO_KEEP,
            )
        )

    # Create ingest config for all stage directories
    stage_dirs = [config.bucket_stage_dir for config in create_ingest_files_configs]

    ingest_config = schemas.component_configs.IngestFilesConfig(
        project_id=settings.GCP_PROJECT_ID,
        nexus_backend_api_url=settings.SITE_URL,
        bigquery_dataset=omics_dataset.name,
        bucket_name=settings.BUCKET_NAME_PRIVATE,
        bucket_stage_dirs=stage_dirs,
        num_workers=settings.INGEST_NUM_WORKERS,
    )

    # Save configs to GCS
    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"

    # Save create_ingest_files configs
    create_ingest_configs_paths = utils.workflows_configs.dump_configs_to_bucket(
        create_ingest_files_configs, configs_stage_dir
    )

    ingest_config_paths = utils.workflows_configs.dump_configs_to_bucket([ingest_config], configs_stage_dir)

    return submit_pipeline(
        pipeline_component=ingest_data_to_bigquery_pipeline,
        display_name=f"Nexus Ingest Data - {omics_dataset.name}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={
            "create_ingest_files_configs": create_ingest_configs_paths,
            "ingest_task_config": ingest_config_paths[0],
        },
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
        labels={"application": settings.GCP_APPLICATION_BILLING_LABEL, "method": "extract"},
    )


def submit_validation_pipeline(
    adata_gcs_paths: list[str],
    validation_report_id: int,
    validation_methods: list[str],
) -> str:
    """
    Submit a Kubeflow pipeline for validating AnnData files.

    Create validation configs in batches (each with max 100 files), save them to GCS,
    and submit the validation pipeline to process them in parallel.

    :param adata_gcs_paths: List of GCS paths to AnnData files to validate
    :param validation_report_id: ID of the validation report to update with results
    :param validation_methods: List of validation method names to apply to each file

    :raise IOError: If there's an error writing config to GCS
    :raise google.api_core.exceptions.GoogleAPIError: If pipeline submission fails

    :return: URL to the Vertex AI Pipeline dashboard for the submitted job
    """
    from cellarium.nexus.workflows.kubeflow.pipelines import validate_anndata_pipeline  # noqa
    from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline  # noqa

    # Create validation configs in batches
    validation_configs = []
    max_files_per_batch = settings.MAX_ADATA_FILES_PER_VALIDATION_BATCH

    # Split the adata_gcs_paths into batches of max_files_per_batch
    for i in range(0, len(adata_gcs_paths), max_files_per_batch):
        batch_paths = adata_gcs_paths[i : i + max_files_per_batch]

        validation_configs.append(
            schemas.component_configs.ValidationConfig(
                nexus_backend_api_url=settings.SITE_URL,
                validation_report_id=validation_report_id,
                adata_gcs_paths=batch_paths,
                validation_methods=validation_methods,
                max_bytes_valid_per_file=settings.INGEST_INPUT_FILE_MAX_SIZE,
            )
        )

    # Save configs to GCS
    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"
    config_paths = utils.workflows_configs.dump_configs_to_bucket(
        configs=validation_configs, bucket_path=configs_stage_dir
    )

    # Submit pipeline and return the pipeline URL
    return submit_pipeline(
        pipeline_component=validate_anndata_pipeline,
        display_name=f"Nexus Validate AnnData - Report {validation_report_id}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={"validation_configs": config_paths},
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
        labels={"application": settings.GCP_APPLICATION_BILLING_LABEL},
    )
