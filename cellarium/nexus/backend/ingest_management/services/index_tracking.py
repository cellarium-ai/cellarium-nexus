from typing import Type

from django.contrib.contenttypes.models import ContentType
from django.db import transaction
from django.db.models import Model
from cellarium.nexus.backend.ingest_management import models


def reserve_indexes(model_class: Type[Model], batch_size: int) -> tuple[int, int]:
    """
    Reserves a batch of indexes atomically for the given model.

    :param model_class: The Django model class to track indexes for.
    :param batch_size: The number of indexes to reserve.

    :return: Tuple (start_index, end_index)
    """
    with transaction.atomic():
        content_type = ContentType.objects.get_for_model(model=model_class)
        tracking, created = models.IndexTracking.objects.select_for_update().get_or_create(content_type=content_type)

        start_index = tracking.largest_index + 1
        end_index = start_index + batch_size - 1

        tracking.largest_index = end_index
        tracking.save()

    return start_index, end_index
