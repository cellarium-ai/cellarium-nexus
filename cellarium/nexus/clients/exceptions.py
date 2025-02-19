class ClientBaseError(Exception):
    pass


class HTTPError(ClientBaseError):
    pass


class HTTPError5XX(HTTPError):
    pass


class HTTPError404(HTTPError):
    pass


class HTTPError403(HTTPError):
    pass


class HTTPError401(HTTPError):
    pass


class HTTPClientError(HTTPError):
    pass
