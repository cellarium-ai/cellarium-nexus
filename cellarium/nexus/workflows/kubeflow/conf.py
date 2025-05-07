BASE_IMAGE: str | None = None
SERVICE_ACCOUNT: str | None = None
LABELS: dict[str, str] = {}


def set_configs(base_image: str, service_account: str, labels: dict[str, str] | None = None):
    """
    Set the configuration values for the global variables.

    This function updates the variables to the provided inputs, allowing for the configuration values to be used
    throughout the application where needed.

    :param base_image: A string representing the base image to be used.
    :param service_account: A string representing the service account to use for authentication.
    :param labels: A dictionary of labels to be applied to the vertex AI components.
    """
    global BASE_IMAGE, SERVICE_ACCOUNT, LABELS
    BASE_IMAGE = base_image
    SERVICE_ACCOUNT = service_account
    LABELS = labels or {}
