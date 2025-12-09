def get_partition_slice(
    total_items: int,
    partition_index: int,
    num_partitions: int,
) -> tuple[int, int]:
    """
    Calculate [start, end) slice for this partition.
    Distributes the remainder over the first `remainder` partitions.

    :param total_items: Total number of items in the sequence
    :param partition_index: Index for which slice has to be calculated
    :param num_partitions: Total number of partitions

    :return: Tuple with start end slice indexes for current partition
    """
    if not (0 <= partition_index < num_partitions):
        raise ValueError(f"partition_index {partition_index} must be in [0, {num_partitions})")
    if num_partitions <= 0:
        raise ValueError(f"num_partitions must be positive, got {num_partitions}")
    if total_items < 0:
        raise ValueError(f"total_items must be non-negative, got {total_items}")

    base = total_items // num_partitions
    rem = total_items % num_partitions

    if partition_index < rem:
        # Partitions [0, rem) get base + 1 items
        start = partition_index * (base + 1)
        end = start + (base + 1)
    else:
        # Remaining partitions get base items
        start = rem * (base + 1) + (partition_index - rem) * base
        end = start + base

    return start, min(end, total_items)
