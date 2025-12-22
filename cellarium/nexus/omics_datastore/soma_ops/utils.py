def get_block_slice(
    total_items: int,
    partition_index: int,
    block_size: int,
) -> tuple[int, int]:
    """
    Calculate [start, end) slice for a fixed-size block.
    Each partition gets at most `block_size` items.
    All full blocks are size `block_size`, the last block gets the remainder.

    :param total_items: Total number of items in the sequence
    :param partition_index: Zero-based partition index
    :param block_size: Target number of items per partition (e.g. 500)

    :return: (start, end) slice indexes for this partition
    """
    if block_size <= 0:
        raise ValueError(f"block_size must be positive, got {block_size}")
    if total_items < 0:
        raise ValueError(f"total_items must be non-negative, got {total_items}")
    if partition_index < 0:
        raise ValueError(f"partition_index must be non-negative, got {partition_index}")

    start = partition_index * block_size
    if start >= total_items:
        # No work for this partition
        return total_items, total_items

    end = min(start + block_size, total_items)
    return start, end
