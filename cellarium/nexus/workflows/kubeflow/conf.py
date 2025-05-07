BASE_IMAGE: str | None = None
SERVICE_ACCOUNT: str | None = None


def set_configs(base_image: str, service_account: str):
    """
    Sets the configuration values for the global variables `BASE_IMAGE` and
    `SERVICE_ACCOUNT`. This function updates these variables to the provided inputs,
    allowing for the configuration values to be used throughout the application
    where needed.

    :param base_image: A string representing the base image to be used.
    :param service_account: A string representing the service account to use for
        authentication.
    :return: None
    """
    global BASE_IMAGE, SERVICE_ACCOUNT
    BASE_IMAGE = base_image
    SERVICE_ACCOUNT = service_account
