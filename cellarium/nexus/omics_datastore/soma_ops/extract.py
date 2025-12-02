"""
SOMA data extraction utilities.

This module provides functions to extract data from SOMA experiments into AnnData files.
"""

import logging
import multiprocessing
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Literal

import anndata
import numpy as np
import tiledbsoma
from anndata._core.aligned_df import ImplicitModificationWarning
from anndata.experimental import AnnCollection
from scipy.sparse import coo_matrix
from tqdm import tqdm

from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaExtractError
from cellarium.nexus.shared.schemas.omics_datastore import SomaExtractPlan, SomaJoinIdRange

logger = logging.getLogger(__name__)


def child_init(log_level: str, verbose: bool = True) -> None:
    """
    Configure logging for child processes.

    :param log_level: Logging level to set
    :param verbose: If False, suppress INFO level logging in child processes
    """
    # Suppress ImplicitModificationWarning from anndata
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

    # If verbose is False, set level to WARNING to suppress INFO logs
    effective_level = log_level if verbose else "WARNING"

    logging.basicConfig(
        level=effective_level,
        format="[%(asctime)s] [%(processName)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def extract_range_to_anndata(
    *,
    experiment_uri: str,
    value_filter: str,
    joinid_range: SomaJoinIdRange,
    output_path: Path,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    var_joinids: list[int] | None = None,
    x_layer: str = "X",
    output_format: Literal["zarr", "h5ad"] = "zarr",
) -> None:
    """
    Extract a single soma_joinid range to an AnnData file.

    Read obs, var, and X data from a SOMA experiment for a specific joinid range
    and write to an AnnData file in the specified format.

    :param experiment_uri: URI of the SOMA experiment
    :param value_filter: SOMA obs value_filter expression
    :param joinid_range: Inclusive soma_joinid range to extract
    :param output_path: Local path to save AnnData file
    :param obs_columns: Optional obs columns to include
    :param var_columns: Optional var columns to include
    :param var_joinids: Optional list of var soma_joinids to filter features by
    :param x_layer: Name of the SOMA X layer to read counts from
    :param output_format: Output format - "zarr" or "h5ad"

    :raise SomaExtractError: If SOMA reads fail
    :raise IOError: If file operations fail
    :raise ValueError: If output_format is invalid
    """
    try:
        logger.info(f"Extracting joinid range [{joinid_range.start}, {joinid_range.end}] to {output_path}")

        # Open SOMA experiment
        logger.info(f"Opening SOMA experiment: {experiment_uri}")

        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            # Read obs data
            logger.info(f"Reading obs data for range [{joinid_range.start}, {joinid_range.end}]")

            obs_query = exp.obs.read(
                coords=(slice(joinid_range.start, joinid_range.end),),
                value_filter=value_filter if value_filter else None,
            )

            logger.info("Converting obs query to pandas DataFrame")

            # Convert to pandas DataFrame
            obs_df = obs_query.concat().to_pandas()

            logger.info(f"Found {len(obs_df)} cells in obs")

            if obs_df.empty:
                logger.warning(f"No cells found in joinid range [{joinid_range.start}, {joinid_range.end}]")
                # Create empty AnnData
                adata = anndata.AnnData()
                adata.write_h5ad(output_path, compression="gzip")
                return

            # Filter obs columns if specified
            if obs_columns is not None:
                # Always keep soma_joinid
                columns_to_keep = ["soma_joinid"] + [
                    col for col in obs_columns if col in obs_df.columns and col != "soma_joinid"
                ]
                obs_df = obs_df[columns_to_keep]

            # Set index to soma_joinid
            obs_df = obs_df.set_index("soma_joinid")

            # Read var data
            logger.info("Reading var data")

            rna = exp.ms["RNA"]
            if var_joinids is not None:
                # Filter var by provided joinids, preserving order
                var_joinids_array = np.array(var_joinids, dtype=np.int64)
                var_df = rna.var.read(coords=(var_joinids_array,)).concat().to_pandas()
                var_df = var_df.set_index("soma_joinid")
                # Reindex to preserve the order of var_joinids
                var_df = var_df.reindex(var_joinids_array)
                logger.info(f"Filtered to {len(var_df)} features by var_joinids")
            else:
                var_df = rna.var.read().concat().to_pandas()
                var_df = var_df.set_index("soma_joinid")
                logger.info(f"Found {len(var_df)} features in var")

            var_joinids_array = var_df.index.to_numpy(dtype=np.int64)

            # Filter var columns if specified
            if var_columns is not None:
                # Keep only requested columns (soma_joinid is already the index)
                columns_to_keep = [col for col in var_columns if col in var_df.columns]
                var_df = var_df[columns_to_keep]

            # Read X data
            logger.info(f"Accessing X layer: {x_layer}")

            x_data = rna.X[x_layer]

            # Get the joinids we actually have in obs (ensure int64 for SOMA)
            obs_joinids = np.asarray(obs_df.index.to_numpy(), dtype=np.int64)

            logger.info(f"Reading X data for {len(obs_joinids)} cells and {len(var_joinids_array)} features")

            # Read X data for specific obs and var joinids
            x_query = x_data.read(
                coords=(obs_joinids, var_joinids_array),
            )

            logger.info("Converting X query to COO format")

            # Correct pattern for SparseNDArrayRead → Arrow SparseCOOTensor → SciPy
            x_tensor = x_query.coos().concat()
            x_coo = x_tensor.to_scipy()

            logger.info(f"X data shape: {x_coo.shape}, nnz: {x_coo.nnz}")

            # Convert to coo_matrix (anndata doesn't accept coo_array)
            logger.info("Converting to COO matrix and remapping indices")

            x_coo_matrix = coo_matrix(x_coo)

            # Remap row and column indices from global soma_joinid space to local 0-based indices
            # Use numpy searchsorted for vectorized O(n log n) lookup instead of Python dict O(n)
            obs_joinids_sorted_idx = np.argsort(obs_joinids)
            obs_joinids_sorted = obs_joinids[obs_joinids_sorted_idx]
            row_positions = np.searchsorted(obs_joinids_sorted, x_coo_matrix.row)
            remapped_rows = obs_joinids_sorted_idx[row_positions]

            var_joinids_sorted_idx = np.argsort(var_joinids_array)
            var_joinids_sorted = var_joinids_array[var_joinids_sorted_idx]
            col_positions = np.searchsorted(var_joinids_sorted, x_coo_matrix.col)
            remapped_cols = var_joinids_sorted_idx[col_positions]

            # Create new sparse matrix with remapped coordinates in CSR format
            # (AnnData prefers CSR for storage and can't write COO to HDF5)
            x_coo_remapped = coo_matrix(
                (x_coo_matrix.data, (remapped_rows, remapped_cols)), shape=(len(obs_df), len(var_df))
            )
            x_csr = x_coo_remapped.tocsr()

            # Create AnnData object
            logger.info("Creating AnnData object")

            adata = anndata.AnnData(X=x_csr, obs=obs_df, var=var_df)

            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Write to file
            logger.info(f"Writing to {output_path}")

            if output_format == "zarr":
                adata.write_zarr(output_path)
            elif output_format == "h5ad":
                adata.write_h5ad(output_path, compression="gzip")
            else:
                raise ValueError(f"output_format must be 'zarr' or 'h5ad', got {output_format}")

            logger.info(f"Successfully extracted {len(obs_df)} cells to {output_path}")

    except Exception as e:
        logger.error(f"Failed to extract joinid range [{joinid_range.start}, {joinid_range.end}]: {e}")
        raise SomaExtractError(f"SOMA extraction failed for range [{joinid_range.start}, {joinid_range.end}]") from e


def _extract_range_worker(
    idx: int,
    experiment_uri: str,
    value_filter: str,
    joinid_range: SomaJoinIdRange,
    output_path: Path,
    obs_columns: list[str] | None,
    var_columns: list[str] | None,
    var_joinids: list[int] | None,
    x_layer: str,
    output_format: Literal["zarr", "h5ad"],
) -> tuple[int, str]:
    """
    Worker function for parallel extraction.

    :param idx: Index of the range being processed
    :param experiment_uri: URI of the SOMA experiment
    :param value_filter: SOMA value filter string for obs
    :param joinid_range: Range of soma_joinids to extract
    :param output_path: Path to write the output file
    :param obs_columns: Optional list of obs columns to include
    :param var_columns: Optional list of var columns to include
    :param var_joinids: Optional list of var soma_joinids to filter features by
    :param x_layer: Name of the X layer to read
    :param output_format: Output format, either "zarr" or "h5ad"

    :raise SomaExtractError: If extraction fails

    :return: Tuple of (range_index, output_path_str)
    """
    extract_range_to_anndata(
        experiment_uri=experiment_uri,
        value_filter=value_filter,
        joinid_range=joinid_range,
        output_path=output_path,
        obs_columns=obs_columns,
        var_columns=var_columns,
        var_joinids=var_joinids,
        x_layer=x_layer,
        output_format=output_format,
    )

    return idx, str(output_path)


def extract_ranges(
    *,
    plan: SomaExtractPlan,
    output_dir: Path,
    output_format: Literal["zarr", "h5ad"] = "zarr",
    max_workers: int | None = None,
    verbose: bool = True,
) -> None:
    """
    Extract multiple soma_joinid ranges in parallel.

    Process all joinid ranges from the extract plan and write separate AnnData files
    for each range in the specified format. All data specification (obs_columns,
    var_columns, var_joinids, x_layer) is taken from the plan.

    :param plan: SOMA extract plan with all data specification
    :param output_dir: Local directory to save AnnData files
    :param output_format: Output format - "zarr" or "h5ad" (default: "zarr")
    :param max_workers: Maximum number of parallel workers
    :param verbose: If False, suppress INFO level logging in parallel workers

    :raise SomaExtractError: If SOMA reads fail
    :raise IOError: If file operations fail
    :raise ValueError: If output_format is invalid
    """
    # Suppress ImplicitModificationWarning from anndata
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Extracting {len(plan.joinid_ranges)} ranges using {max_workers} workers (spawn mode)")

    # Execute in parallel
    completed = 0
    failed = []

    # Use spawn method to avoid fork-related issues with network connections (TileDB/SOMA)
    mp_context = multiprocessing.get_context("spawn")

    # Determine file extension based on output format
    file_ext = ".zarr" if output_format == "zarr" else ".h5ad"

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp_context,
        initializer=child_init,
        initargs=(logging.getLevelName(logger.level), verbose),
    ) as executor:
        # Submit all extraction jobs
        futures = {}
        for idx, joinid_range in enumerate(plan.joinid_ranges):
            output_path = output_dir / f"range_{idx:06d}{file_ext}"
            future = executor.submit(
                _extract_range_worker,
                idx,
                plan.experiment_uri,
                plan.value_filter,
                joinid_range,
                output_path,
                plan.obs_columns,
                plan.var_columns,
                plan.var_joinids,
                plan.x_layer,
                output_format,
            )
            futures[future] = idx

        # Wait for completion
        for future in tqdm(as_completed(futures), total=len(plan.joinid_ranges), desc="Extracting ranges"):
            try:
                idx, output_path = future.result()
                completed += 1
            except Exception as e:
                idx = futures[future]
                failed.append((idx, str(e)))
                logger.error(f"Range {idx} failed: {e}")

    if failed:
        error_msg = f"Failed to extract {len(failed)} ranges: {failed}"
        logger.error(error_msg)
        raise SomaExtractError(error_msg)

    logger.info(f"Successfully extracted all {len(plan.joinid_ranges)} ranges to {output_dir}")


def _write_shuffled_chunk(
    chunk_idx: int,
    chunk_indices: np.ndarray,
    input_files: list[Path],
    input_format: Literal["zarr", "h5ad"],
    output_dir: Path,
    output_format: Literal["zarr", "h5ad"],
) -> tuple[int, str]:
    """
    Worker function to write a single shuffled chunk.

    Each worker creates its own AnnCollection from all input files to enable
    complete shuffling across all ranges while avoiding pickling issues.

    :param chunk_idx: Index of the chunk being written
    :param chunk_indices: Array of cell indices for this chunk
    :param input_files: List of all input file paths
    :param input_format: Input format
    :param output_dir: Output directory
    :param output_format: Output format

    :return: Tuple of (chunk_idx, output_path_str)
    """
    logger.info(f"Processing chunk {chunk_idx} ({len(chunk_indices)} cells)...")

    # Each worker creates its own AnnCollection with ALL files (backed mode)
    if input_format == "zarr":
        adatas = [anndata.read_zarr(str(f)) for f in input_files]
    else:  # h5ad
        adatas = [anndata.read_h5ad(str(f), backed="r") for f in input_files]

    ann_collection = AnnCollection(
        adatas,
        join_obs="inner",
        join_vars="inner",
    )

    # Extract shuffled cells (may come from any/all ranges)
    chunk_adata = ann_collection[chunk_indices].to_adata()

    # Write output
    if output_format == "zarr":
        output_path = output_dir / f"chunk_{chunk_idx:06d}.zarr"
        chunk_adata.write_zarr(output_path)
    else:  # h5ad
        output_path = output_dir / f"chunk_{chunk_idx:06d}.h5ad"
        chunk_adata.write_h5ad(output_path, compression="gzip")

    logger.info(f"Successfully wrote chunk {chunk_idx} with {len(chunk_indices)} cells to {output_path}")

    return chunk_idx, str(output_path)


def shuffle_extracted_chunks(
    *,
    input_dir: Path,
    output_dir: Path,
    chunk_size: int,
    input_format: Literal["zarr", "h5ad"] = "zarr",
    output_format: Literal["zarr", "h5ad"] = "h5ad",
    max_workers: int | None = None,
    verbose: bool = True,
) -> None:
    """
    Shuffle cells across extracted AnnData chunks.

    Read contiguous range files from input directory, redistribute cells
    randomly across new chunks of the specified size, and write shuffled
    output in the specified format.

    :param input_dir: Directory with contiguous range files
    :param output_dir: Directory to write shuffled chunks
    :param chunk_size: Number of cells per output chunk
    :param input_format: Input format - "zarr" or "h5ad" (default: "zarr")
    :param output_format: Output format - "zarr" or "h5ad" (default: "h5ad")
    :param max_workers: Maximum parallel workers for writing
    :param verbose: If False, suppress INFO level logging in parallel workers

    :raise IOError: If file operations fail
    :raise ValueError: If chunk_size is not positive or formats are invalid
    """
    # Suppress ImplicitModificationWarning from anndata
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    if input_format not in ("zarr", "h5ad"):
        raise ValueError(f"input_format must be 'zarr' or 'h5ad', got {input_format}")

    if output_format not in ("zarr", "h5ad"):
        raise ValueError(f"output_format must be 'zarr' or 'h5ad', got {output_format}")

    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    logger.info(f"Shuffling cells from {input_dir} ({input_format}) to {output_dir} ({output_format})")

    # Find all input files based on input format
    input_pattern = f"range_*.{input_format}" if input_format == "h5ad" else "range_*.zarr"
    input_files = sorted(input_dir.glob(input_pattern))
    if not input_files:
        raise ValueError(f"No range files found in {input_dir} with pattern {input_pattern}")

    logger.info(f"Found {len(input_files)} input files")

    # Count total cells by reading metadata from files
    logger.info("Computing total cell count...")
    if input_format == "zarr":
        total_cells = sum(anndata.read_zarr(str(f)).n_obs for f in input_files)
    else:  # h5ad
        total_cells = sum(anndata.read_h5ad(str(f), backed="r").n_obs for f in input_files)

    logger.info(f"Total cells to shuffle: {total_cells}")

    # Create random permutation of all cell indices
    logger.info("Computing random permutation...")
    cell_permutation = np.random.permutation(total_cells)

    # Compute output chunks
    num_chunks = (total_cells + chunk_size - 1) // chunk_size
    logger.info(f"Creating {num_chunks} output chunks")

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process each output chunk in parallel
    logger.info(f"Writing shuffled chunks using {max_workers} workers (spawn mode)...")

    # Use spawn method for multiprocessing
    mp_context = multiprocessing.get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp_context,
        initializer=child_init,
        initargs=(logging.getLevelName(logger.level), verbose),
    ) as executor:
        futures = {}
        for chunk_idx in range(num_chunks):
            start_idx = chunk_idx * chunk_size
            end_idx = min(start_idx + chunk_size, total_cells)
            chunk_indices = cell_permutation[start_idx:end_idx]

            # Submit chunk processing job
            # Each worker will create its own AnnCollection with all files
            future = executor.submit(
                _write_shuffled_chunk,
                chunk_idx,
                chunk_indices,
                input_files,
                input_format,
                output_dir,
                output_format,
            )
            futures[future] = chunk_idx

        # Wait for completion and track progress
        completed = 0
        for future in tqdm(as_completed(futures), total=num_chunks, desc="Shuffling chunks"):
            try:
                chunk_idx, output_path = future.result()
                completed += 1
            except Exception as e:
                chunk_idx = futures[future]
                logger.error(f"Failed to write chunk {chunk_idx}: {e}")
                raise

    logger.info(f"Successfully shuffled {total_cells} cells into {num_chunks} chunks")
