import functools
import os
import tempfile
import typing as t

from google.cloud import aiplatform
from google_cloud_pipeline_components.v1.custom_job import create_custom_training_job_from_component
from kfp import compiler, dsl
from kfp.components import BaseComponent

from cellarium.nexus.workflows.kubeflow.utils import constants


def dsl_component_job(
    base_image: str,
    display_name: str = "",
    machine_type: str = "n1-standard-4",
    replica_count: int = 1,
    **custom_job_params,
) -> t.Callable[[t.Callable], t.Callable]:
    def decorator(func: t.Callable) -> t.Callable:
        component = dsl.component(base_image=base_image)(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            job = create_custom_training_job_from_component(
                component_spec=component,
                display_name=display_name,
                machine_type=machine_type,
                replica_count=replica_count,
                **custom_job_params,
            )(*args, **kwargs)

            return job

        return wrapper

    return decorator


def submit_pipeline(
    pipeline_component: BaseComponent | t.Callable,
    display_name: str,
    gcp_project: str,
    pipeline_kwargs: dict[str, t.Any],
    pipeline_location: str = constants.DEFAULT_PIPELINE_LOCATION,
    service_account: str | None = None,
    pipeline_root_path: str | None = None,
    labels: dict[str, str] | None = None,
) -> str:
    """
    Create and run a pipeline on Vertex AI Pipelines. Use a temporary file to compile the pipeline config,
    then run the pipeline job and delete the temporary file.

    :param pipeline_component: Pipeline function, must be a function wrapped :func:`kfp.dsl.pipeline` decorator.
    :param display_name: A name displayed in the Vertex AI Pipelines UI.
    :param gcp_project: GCP Project where the pipeline should be submitted to.
    :param pipeline_kwargs: Keyword arguments to pass to the pipeline function.
    :param pipeline_location: Datacenter location of Google Cloud Platform to run the pipeline job.
    :param service_account: Service account email to use for the pipeline execution. If None, the default compute
        service account will be used.
    :param pipeline_root_path: GCS path to use as the pipeline root directory. If None, the default
        Vertex AI Pipelines location will be used.
    :param labels: Labels to apply to the pipeline job.

    :raise: google.api_core.exceptions.GoogleAPIError: If pipeline submission fails.

    :return: URL link to the Vertex AI Pipeline dashboard for the submitted job
    """
    temp_file = tempfile.NamedTemporaryFile(suffix=".yaml")
    os.environ["GRPC_DNS_RESOLVER"] = "native"

    aiplatform.init(location=pipeline_location, project=gcp_project)

    compiler.Compiler().compile(pipeline_func=pipeline_component, package_path=temp_file.name)

    job = aiplatform.PipelineJob(
        display_name=display_name,
        template_path=temp_file.name,
        parameter_values=pipeline_kwargs,
        pipeline_root=pipeline_root_path,
        labels=labels,
    )

    job.submit(service_account=service_account)
    temp_file.close()

    # Generate the Vertex AI Pipeline dashboard URL
    return (
        f"https://console.cloud.google.com/vertex-ai/locations/{pipeline_location}/"
        f"pipelines/runs/{job.name}?project={gcp_project}"
    )
