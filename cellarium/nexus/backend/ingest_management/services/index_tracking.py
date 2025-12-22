from django.db import transaction

from cellarium.nexus.backend.ingest_management import models as ingest_models


def reserve_indexes(
    *, omics_dataset: "ingest_models.OmicsDataset", resource_key: str, batch_size: int
) -> tuple[int, int]:
    """
    Reserve a batch of indexes atomically for the given dataset and resource key.

    :param omics_dataset: Omics dataset instance to scope the index space

    :param resource_key: Logical resource key identifying the index namespace (e.g., "cell_info")

    :param batch_size: Number of indexes to reserve in this batch

    :return: Tuple with start and end indexes (start_index, end_index)
    """
    with transaction.atomic():
        tracking, _ = ingest_models.IndexTracking.objects.select_for_update().get_or_create(
            omics_dataset=omics_dataset,
            resource_key=resource_key,
        )

        start_index = tracking.largest_index + 1
        end_index = start_index + batch_size - 1

        tracking.largest_index = end_index
        tracking.save()

    return start_index, end_index
