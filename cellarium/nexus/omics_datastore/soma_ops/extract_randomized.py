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
from scipy.sparse import coo_matrix
from tenacity import before_log, retry, stop_after_attempt, wait_exponential
from tqdm import tqdm

from cellarium.nexus.omics_datastore.soma_ops import utils
from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaExtractError
from cellarium.nexus.shared.schemas.omics_datastore import IdContiguousRange, RandomizedCurriculumMetadata

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


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    before=before_log(logger, logging.INFO),
    reraise=True,
)
def extract_range_to_anndata(
    *,
    experiment_uri: str,
    value_filter: str,
    joinid_range: IdContiguousRange,
    output_path: Path,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    x_layer: str = "X",
    output_format: Literal["zarr", "h5ad"] = "zarr",
) -> None:
    """
    Extract a single soma_joinid range to an AnnData file.

    Read obs, var, and X data from a SOMA experiment for a specific joinid range
    and write to an AnnData file in the specified format. Feature filtering is
    deferred to the shuffle stage for better SOMA query performance.

    :param experiment_uri: URI of the SOMA experiment
    :param value_filter: SOMA obs value_filter expression
    :param joinid_range: Inclusive soma_joinid range to extract
    :param output_path: Local path to save AnnData file
    :param obs_columns: Optional obs columns to include
    :param var_columns: Optional var columns to include
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

            # Read var data (all features - filtering deferred to shuffle stage)
            logger.info("Reading var data")

            rna = exp.ms["RNA"]
            var_df = rna.var.read().concat().to_pandas()
            var_df = var_df.set_index("soma_joinid")
            logger.info(f"Found {len(var_df)} features in var")

            var_joinids_array = var_df.index.to_numpy(dtype=np.int64)

            # Filter var columns if specified
            if var_columns is not None:
                # Keep only requested columns (soma_joinid is already the index, skip it)
                columns_to_keep = [col for col in var_columns if col in var_df.columns and col != "soma_joinid"]
                missing_cols = [col for col in var_columns if col not in var_df.columns and col != "soma_joinid"]
                if missing_cols:
                    logger.warning(
                        f"Requested var columns not found: {missing_cols}. "
                        f"Available columns: {var_df.columns.tolist()}"
                    )
                if columns_to_keep:
                    var_df = var_df[columns_to_keep]
                    logger.info(f"Keeping var columns: {columns_to_keep}")
                else:
                    logger.info("No matching var columns found, keeping all columns")

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
    joinid_range: IdContiguousRange,
    output_path: Path,
    obs_columns: list[str] | None,
    var_columns: list[str] | None,
    x_layer: str,
    output_format: Literal["zarr", "h5ad"],
) -> tuple[int, str]:
    """
    Worker function for parallel extraction.

    :param idx: Index of the range being processed.
    :param experiment_uri: URI of the SOMA experiment.
    :param value_filter: SOMA value filter string for obs.
    :param joinid_range: Range of soma_joinids to extract.
    :param output_path: Path to write the output file.
    :param obs_columns: Optional list of obs columns to include.
    :param var_columns: Optional list of var columns to include.
    :param x_layer: Name of the X layer to read.
    :param output_format: Output format, either "zarr" or "h5ad".

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
        x_layer=x_layer,
        output_format=output_format,
    )

    return idx, str(output_path)


def extract_ranges(
    *,
    curriculum_metadata: RandomizedCurriculumMetadata,
    output_dir: Path,
    partition_index: int = 0,
    max_ranges_per_partition: int | None = None,
    output_format: Literal["zarr", "h5ad"] = "zarr",
    max_workers: int | None = None,
    verbose: bool = True,
) -> None:
    """
    Extract multiple soma_joinid ranges in parallel.

    Process all joinid ranges from the curriculum metadata and write separate AnnData files
    for each range in the specified format. All data specification (obs_columns,
    var_columns, var_joinids, x_layer) is taken from the curriculum metadata.

    :param curriculum_metadata: SOMA extract metadata with all data specification.
    :param output_dir: Local directory to save AnnData files.
    :param partition_index: Index used for slicing ranges and output chunk indexes.
        Needed for distributing extracting over multiple distributed VMs. Default is 0 (single VM execution).
    :param max_ranges_per_partition: Partition block size. Default is None, this means it will use all ranges
    :param output_format: Output format - "zarr" or "h5ad" (default: "zarr").
    :param max_workers: Maximum number of parallel workers.
    :param verbose: If False, suppress INFO level logging in parallel workers.

    :raise SomaExtractError: If SOMA reads fail
    :raise IOError: If file operations fail
    :raise ValueError: If output_format is invalid
    """
    # Suppress ImplicitModificationWarning from anndata
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

    if max_ranges_per_partition is None:
        max_ranges_per_partition = curriculum_metadata.num_ranges

    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Calculating current extract ranges and output ids
    range_slice_start, range_slice_end = utils.get_block_slice(
        total_items=curriculum_metadata.num_ranges,
        partition_index=partition_index,
        block_size=max_ranges_per_partition,
    )

    ranges_to_process = curriculum_metadata.id_ranges[range_slice_start:range_slice_end]

    logger.info(f"Extracting {len(ranges_to_process)} ranges using {max_workers} workers (spawn mode)")

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
        for idx, joinid_range in enumerate(ranges_to_process):
            output_path = output_dir / f"range_{idx:06d}{file_ext}"
            future = executor.submit(
                _extract_range_worker,
                idx,
                curriculum_metadata.experiment_uri,
                curriculum_metadata.value_filter,
                joinid_range,
                output_path,
                curriculum_metadata.obs_columns,
                curriculum_metadata.var_columns,
                curriculum_metadata.x_layer,
                output_format,
            )
            futures[future] = idx

        # Wait for completion
        for future in tqdm(as_completed(futures), total=len(ranges_to_process), desc="Extracting ranges"):
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

    logger.info(f"Successfully extracted all {len(ranges_to_process)} ranges to {output_dir}")


def consolidate_zarr_extracts(
    *,
    input_dir: Path,
    output_path: Path,
    batch_size: int = 50,
) -> int:
    """
    Consolidate multiple Zarr range files into a single Zarr file.

    Read all range_*.zarr files from input directory and concatenate them
    into a single consolidated Zarr file for efficient random access during shuffling.

    :param input_dir: Directory containing range_*.zarr files
    :param output_path: Path to write the consolidated Zarr file
    :param batch_size: Number of files to concatenate per batch to manage memory

    :raise ValueError: If no input files found
    :raise IOError: If file operations fail

    :return: Total number of cells in the consolidated file
    """
    zarr_files = sorted(input_dir.glob("range_*.zarr"))

    if not zarr_files:
        raise ValueError(f"No range_*.zarr files found in {input_dir}")

    logger.info(f"Consolidating {len(zarr_files)} Zarr files into {output_path}")

    # Concatenate in batches to manage memory
    consolidated = None

    for i in range(0, len(zarr_files), batch_size):
        batch_files = zarr_files[i : i + batch_size]
        batch_adatas = [anndata.read_zarr(str(f)) for f in batch_files]
        batch_concat = anndata.concat(batch_adatas, join="outer")

        if consolidated is None:
            consolidated = batch_concat
        else:
            consolidated = anndata.concat([consolidated, batch_concat], join="outer")

        logger.info(f"Processed {min(i + batch_size, len(zarr_files))}/{len(zarr_files)} files")

        # Free memory
        del batch_adatas
        del batch_concat

    # Write consolidated file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    consolidated.write_zarr(output_path)

    total_cells = consolidated.n_obs
    logger.info(f"Consolidated {total_cells} cells into {output_path}")

    return total_cells


def _write_shuffle_chunk(
    consolidated: anndata.AnnData,
    chunk_indices: np.ndarray,
    output_path: Path,
    var_joinids: list[int] | None,
) -> str:
    """
    Write a single shuffled chunk to h5ad file.

    :param consolidated: Consolidated AnnData object (shared via fork copy-on-write)
    :param chunk_indices: Array of cell indices for this chunk
    :param output_path: Path to write the output file
    :param var_joinids: Optional list of var soma_joinids to filter features by

    :raise IOError: If file operations fail

    :return: Output path as string
    """
    chunk_adata = consolidated[chunk_indices].copy()

    # Apply feature filtering if var_joinids provided
    if var_joinids is not None:
        chunk_adata = chunk_adata[:, var_joinids].copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    chunk_adata.write_h5ad(output_path, compression="gzip")

    return str(output_path)


def shuffle_extracted_chunks(
    *,
    curriculum_metadata: RandomizedCurriculumMetadata,
    input_dir: Path,
    output_dir: Path,
    partition_index: int = 0,
    max_output_chunks_per_partition: int | None = None,
    max_workers: int | None = None,
    consolidate_batch_size: int = 50,
) -> None:
    """
    Shuffle cells across extracted AnnData chunks.

    Consolidate Zarr range files into a single file, then redistribute cells
    randomly across output chunks using fork-based multiprocessing for
    memory-efficient parallel writes.

    :param curriculum_metadata: SOMA extract metadata with all data specification
    :param input_dir: Directory with range_*.zarr files from extraction stage
    :param output_dir: Directory to write shuffled h5ad chunks
    :param partition_index: Index used for slicing output chunk indexes.
        Needed for distributing across multiple distributed VMs. Default is 0.
    :param max_output_chunks_per_partition: Partition block size. Default is None (all chunks)
    :param max_workers: Maximum parallel workers for writing shuffled chunks
    :param consolidate_batch_size: Batch size for consolidation step

    :raise IOError: If file operations fail
    :raise ValueError: If no input files found
    :raise SomaExtractError: If bin count mismatch
    """
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    if max_output_chunks_per_partition is None:
        max_output_chunks_per_partition = curriculum_metadata.num_bins

    # Stage 1: Consolidate Zarr files into single file
    logger.info("Stage 1: Consolidating Zarr files...")
    consolidated_path = input_dir / "_consolidated.zarr"
    total_cells = consolidate_zarr_extracts(
        input_dir=input_dir,
        output_path=consolidated_path,
        batch_size=consolidate_batch_size,
    )

    # Stage 2: Load consolidated file and prepare shuffle
    logger.info("Stage 2: Loading consolidated data and preparing shuffle...")
    consolidated = anndata.read_zarr(str(consolidated_path))

    chunk_size = curriculum_metadata.extract_bin_size
    extract_bin_id_slice_start, extract_bin_id_slice_end = utils.get_block_slice(
        total_items=curriculum_metadata.num_bins,
        partition_index=partition_index,
        block_size=max_output_chunks_per_partition,
    )
    extract_bin_indexes = curriculum_metadata.extract_bin_indexes[extract_bin_id_slice_start:extract_bin_id_slice_end]

    _num_bins_computed = (total_cells + chunk_size - 1) // chunk_size
    num_bins = len(extract_bin_indexes)

    if _num_bins_computed != num_bins:
        raise SomaExtractError(
            f"Number of extract bin indexes doesn't accommodate all cells within this worker. "
            f"Required number of extract bins {_num_bins_computed}, while "
            f"number of indexes for naming {num_bins}"
        )

    # Create random permutation
    cell_permutation = np.random.permutation(total_cells)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Stage 3: Write shuffled chunks in parallel using fork (copy-on-write)
    logger.info(f"Stage 3: Writing {num_bins} shuffled chunks using {max_workers} workers (fork mode)...")

    mp_context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=max_workers, mp_context=mp_context) as executor:
        futures = {}
        for bin_idx in range(num_bins):
            start_idx = bin_idx * chunk_size
            end_idx = min(start_idx + chunk_size, total_cells)
            chunk_indices = cell_permutation[start_idx:end_idx]
            extract_bin_idx = extract_bin_indexes[bin_idx]

            output_path = output_dir / f"extract_{extract_bin_idx:06d}.h5ad"

            future = executor.submit(
                _write_shuffle_chunk,
                consolidated,
                chunk_indices,
                output_path,
                curriculum_metadata.var_joinids,
            )
            futures[future] = bin_idx

        for future in tqdm(as_completed(futures), total=num_bins, desc="Shuffling extracts"):
            try:
                future.result()
            except Exception as e:
                bin_idx = futures[future]
                logger.error(f"Failed to write extract {bin_idx}: {e}")
                raise

    logger.info(f"Successfully shuffled {total_cells} cells into {num_bins} extracts")
