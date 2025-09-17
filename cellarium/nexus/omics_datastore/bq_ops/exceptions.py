class OmicsDatastoreException(Exception):
    pass


class DataValidationError(OmicsDatastoreException):
    pass


class DataProcessingError(OmicsDatastoreException):
    pass


class DataIngestError(DataProcessingError):
    pass


class BigQueryIngestError(Exception):
    pass


class BigQuerySchemaError(BigQueryIngestError):
    pass


class BigQueryLoadError(BigQueryIngestError):
    pass


class BigQueryCommitError(BigQueryIngestError):
    pass
