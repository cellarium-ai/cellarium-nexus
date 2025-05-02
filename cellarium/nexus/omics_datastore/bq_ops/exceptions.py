class OmicsDatastoreException(Exception):
    pass


class DataProcessingError(OmicsDatastoreException):
    pass


class DataIngestError(DataProcessingError):
    pass
