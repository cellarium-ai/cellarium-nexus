from cellarium.nexus.backend import cell_management
from cellarium.nexus.backend.ingest_management.services import index_tracking


def test_reserve_indexes_creates_and_increments(default_dataset: cell_management.models.OmicsDataset) -> None:
    """
    Reserve indexes twice and assert the second batch continues from the previous end.
    """
    start, end = index_tracking.reserve_indexes(omics_dataset=default_dataset, resource_key="cell_info", batch_size=5)
    assert (start, end) == (1, 5)

    # Second call should continue from previous end
    start2, end2 = index_tracking.reserve_indexes(omics_dataset=default_dataset, resource_key="cell_info", batch_size=3)
    assert (start2, end2) == (6, 8)
