import datetime
import os
import secrets

import pandas as pd
from django.conf import settings

from cellarium.nexus.backend.cell_management.models import BigQueryDataset
from cellarium.nexus.shared import utils
from cellarium.nexus.workflows.kubeflow.component_configs import IngestTaskConfig
from cellarium.nexus.workflows.kubeflow.pipelines import ingest_data_pipeline
from cellarium.nexus.workflows.kubeflow.utils.job import submit_pipeline


def submit_ingest_pipeline(
    df_ingest_file_info: pd.DataFrame, bigquery_dataset: BigQueryDataset, column_mapping: dict
) -> None:
    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    base_stage_dir = f"{settings.BACKEND_PIPELINE_DIR}/data-ingests/{timestamp}_{secrets.token_hex(12)}"

    # Create list for combined task configs
    task_configs = []

    for i, row in df_ingest_file_info.iterrows():
        # Extract file name without extension from gcs_file_path
        file_path = row["gcs_file_path"]
        file_name = os.path.basename(file_path)
        file_name_without_ext = os.path.splitext(file_name)[0]

        # Create unique stage dir for this file
        stage_dir = f"{base_stage_dir}/{secrets.token_hex(8)}_{file_name_without_ext[:10]}"

        tag = row["tag"] if "tag" in df_ingest_file_info.columns else None

        # Create combined config for this task
        task_configs.append(
            IngestTaskConfig(
                project_id=settings.GCP_PROJECT_ID,
                nexus_backend_api_url=settings.SITE_URL,
                bigquery_dataset=bigquery_dataset.name,
                data_source_path=file_path,
                bucket_name=settings.BUCKET_NAME_PRIVATE,
                ingest_bucket_path=stage_dir,
                tag=tag,
                metadata_columns=column_mapping,
            )
        )

    # Save configs to GCS
    configs_stage_dir = f"gs://{settings.BUCKET_NAME_PRIVATE}/pipeline-configs"

    task_config_paths = utils.workflows_configs.dump_configs_to_bucket(task_configs, configs_stage_dir)

    # Submit pipeline
    submit_pipeline(
        pipeline_component=ingest_data_pipeline,
        display_name=f"Nexus Ingest Data - {bigquery_dataset.name}",
        gcp_project=settings.GCP_PROJECT_ID,
        pipeline_kwargs={
            "ingest_task_configs": task_config_paths,
        },
        service_account=settings.PIPELINE_SERVICE_ACCOUNT,
        pipeline_root_path=settings.PIPELINE_ROOT_PATH,
    )
