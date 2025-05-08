class NexusDataOpsError(Exception):
    pass


class NexusDataOpsIngestError(NexusDataOpsError):
    pass


class NexusDataOpsValidationError(NexusDataOpsError):
    pass
