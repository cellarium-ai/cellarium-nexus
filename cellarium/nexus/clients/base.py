import json
from typing import Any

import requests
from nexus.clients import exceptions, status_codes

JsonObject = dict[str, Any]
JsonData = JsonObject | list[JsonObject]


class BaseAPIHTTPClient:
    """
    Base class for communicating with APIs through HTTP protocol

    :param api_url: URL of the Cellarium Cloud Backend API service
    """

    def __init__(self, api_url: str):
        self.api_url = api_url
        super().__init__()

    def _get_endpoint_url(self, endpoint: str) -> str:
        """
        Configure a specific method endpoint from backend url and endpoint

        :param endpoint: Endpoint string without a leading slash
        :return: Full url with backend domains/subdomains and endpoint joined
        """
        return f"{self.api_url}/{endpoint}"

    def _get_headers(self) -> dict[str, str]:
        """
        Get the headers to include in the request

        :return: Headers dictionary
        """
        headers = {}

        # if self.api_token is not None:
        #     headers = {header_names.AUTHORIZATION: f"Bearer {self.api_token}"}

        return headers

    @staticmethod
    def _raise_response_exception(status_code: int, detail: str) -> None:
        """
        Raise an exception based on the status code returned by the server, including the detail message

        :param status_code: HTTP status code
        :param detail: Detail message returned by the server
        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        message = f"Server returned status code {status_code}, Detail: {detail}"
        if status_code == status_codes.STATUS_401_UNAUTHORIZED:
            raise exceptions.HTTPError401(message)
        elif status_code == status_codes.STATUS_403_FORBIDDEN:
            raise exceptions.HTTPError403(message)
        elif status_code == status_codes.STATUS_NOT_FOUND:
            raise exceptions.HTTPError404(message)

        elif (
            status_codes.STATUS_500_INTERNAL_SERVER_ERROR
            <= status_code
            <= status_codes.STATUS_511_NETWORK_AUTHENTICATION_REQUIRED
        ):
            raise exceptions.HTTPError5XX(message)
        else:
            raise exceptions.HTTPError(message)

    def _validate_requests_response(self, response: requests.Response) -> None:
        """
        Validate requests response and raise an exception if response status code is not 200

        :param response: Response object

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        status_code = response.status_code
        if not (status_codes.STATUS_200_OK <= status_code <= status_codes.STATUS_226_IM_USED):
            # When response status code is not 2XX
            try:
                response_detail = response.json()["detail"]
            except (json.decoder.JSONDecodeError, KeyError):
                response_detail = response.text

            self._raise_response_exception(status_code=status_code, detail=response_detail)

    def _request_json(
        self, method: str, endpoint: str, data: JsonData | None = None
    ) -> dict[str, Any] | list[dict[str, Any]]:
        """
        Make a request to the backend service and return JSON response.

        :param method: HTTP method (GET, POST, PUT, PATCH)
        :param endpoint: Endpoint string without a leading slash
        :param data: Payload for POST, PUT, PATCH requests

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
        url = self._get_endpoint_url(endpoint=endpoint)
        headers = self._get_headers()
        response = requests.request(method=method, url=url, headers=headers, json=data)
        self._validate_requests_response(response=response)
        return response.json()

    def get_json(self, endpoint: str) -> JsonData:
        return self._request_json(method="GET", endpoint=endpoint)

    def post_json(self, endpoint: str, data: JsonData | None = None) -> JsonData:
        return self._request_json(method="POST", endpoint=endpoint, data=data)

    def put_json(self, endpoint: str, data: JsonData | None = None) -> JsonData:
        return self._request_json(method="PUT", endpoint=endpoint, data=data)

    def patch_json(self, endpoint: str, data: JsonData | None = None) -> JsonData:
        return self._request_json(method="PATCH", endpoint=endpoint, data=data)
