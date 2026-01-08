"""
SOMA grouped data extraction utilities.

This module provides functions to _extract data from SOMA experiments into AnnData files
where cells are grouped by specified obs columns.
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

from cellarium.nexus.omics_datastore.soma_ops._extract import utils
from cellarium.nexus.omics_datastore.soma_ops.exceptions import SomaExtractError
from cellarium.nexus.shared.schemas.omics_datastore import GroupedBin, GroupedCurriculumMetadata

logger = logging.getLogger(__name__)


def child_init(log_level: str, verbose: bool = True) -> None:
    """
    Configure logging for child processes.

    :param log_level: Logging level to set
    :param verbose: If False, suppress INFO level logging in child processes
    """
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

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
def extract_grouped_bin_to_anndata(
    *,
    experiment_uri: str,
    grouped_bin: GroupedBin,
    output_path: Path,
    value_filter: str | None = None,
    obs_columns: list[str] | None = None,
    var_columns: list[str] | None = None,
    var_joinids: list[int] | None = None,
    x_layer: str = "X",
    output_format: Literal["zarr", "h5ad"] = "h5ad",
) -> None:
    """
    Extract a single grouped bin to an AnnData file.

    Query by joinid range + combined filter (user filter AND group filter) for efficiency.

    :param experiment_uri: URI of the SOMA experiment
    :param grouped_bin: GroupedBin with group filter and joinid bounds
    :param output_path: Local path to save AnnData file
    :param value_filter: Optional user-provided filter to combine with group filter
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
        logger.info(
            f"Extracting grouped bin [{grouped_bin.joinid_min}, {grouped_bin.joinid_max}] "
            f"group_key={grouped_bin.group_key} to {output_path}"
        )

        with tiledbsoma.open(experiment_uri, mode="r") as exp:
            # Combine user filter with group filter
            if value_filter:
                combined_filter = f"({value_filter}) and ({grouped_bin.group_filter})"
            else:
                combined_filter = grouped_bin.group_filter

            logger.info(f"Reading obs data with filter: {combined_filter}")

            obs_query = exp.obs.read(
                coords=(slice(grouped_bin.joinid_min, grouped_bin.joinid_max),),
                value_filter=combined_filter,
            )

            obs_df = obs_query.concat().to_pandas()

            logger.info(f"Found {len(obs_df)} cells in grouped bin")

            if obs_df.empty:
                logger.warning(f"No cells found in grouped bin {grouped_bin.group_key}")
                adata = anndata.AnnData()
                adata.write_h5ad(output_path, compression="gzip")
                return

            # Filter obs columns if specified
            if obs_columns is not None:
                columns_to_keep = ["soma_joinid"] + [
                    col for col in obs_columns if col in obs_df.columns and col != "soma_joinid"
                ]
                obs_df = obs_df[columns_to_keep]

            obs_df = obs_df.set_index("soma_joinid")

            # Read var data
            logger.info("Reading var data")

            rna = exp.ms["RNA"]
            var_df = rna.var.read().concat().to_pandas()
            var_df = var_df.set_index("soma_joinid")
            logger.info(f"Found {len(var_df)} features in var")

            # Apply var filtering if specified
            if var_joinids is not None:
                var_df = var_df.loc[var_joinids]
                logger.info(f"Filtered to {len(var_df)} features")

            var_joinids_array = var_df.index.to_numpy(dtype=np.int64)

            # Filter var columns if specified
            if var_columns is not None:
                columns_to_keep = [col for col in var_columns if col in var_df.columns and col != "soma_joinid"]
                missing_cols = [col for col in var_columns if col not in var_df.columns and col != "soma_joinid"]
                if missing_cols:
                    logger.warning(
                        f"Requested var columns not found: {missing_cols}. "
                        f"Available columns: {var_df.columns.tolist()}"
                    )
                if columns_to_keep:
                    var_df = var_df[columns_to_keep]

            # Read X data
            logger.info(f"Accessing X layer: {x_layer}")

            x_data = rna.X[x_layer]

            obs_joinids = np.asarray(obs_df.index.to_numpy(), dtype=np.int64)

            logger.info(f"Reading X data for {len(obs_joinids)} cells and {len(var_joinids_array)} features")

            x_query = x_data.read(
                coords=(obs_joinids, var_joinids_array),
            )

            logger.info("Converting X query to COO format")

            x_tensor = x_query.coos().concat()
            x_coo = x_tensor.to_scipy()

            logger.info(f"X data shape: {x_coo.shape}, nnz: {x_coo.nnz}")

            # Convert to coo_matrix and remap indices
            logger.info("Converting to COO matrix and remapping indices")

            x_coo_matrix = coo_matrix(x_coo)

            # Remap row indices
            obs_joinids_sorted_idx = np.argsort(obs_joinids)
            obs_joinids_sorted = obs_joinids[obs_joinids_sorted_idx]
            row_positions = np.searchsorted(obs_joinids_sorted, x_coo_matrix.row)
            remapped_rows = obs_joinids_sorted_idx[row_positions]

            # Remap column indices
            var_joinids_sorted_idx = np.argsort(var_joinids_array)
            var_joinids_sorted = var_joinids_array[var_joinids_sorted_idx]
            col_positions = np.searchsorted(var_joinids_sorted, x_coo_matrix.col)
            remapped_cols = var_joinids_sorted_idx[col_positions]

            # Create sparse matrix with remapped coordinates
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
        logger.error(f"Failed to _extract grouped bin {grouped_bin.group_key}: {e}")
        raise SomaExtractError(f"SOMA extraction failed for grouped bin {grouped_bin.group_key}") from e


def _extract_grouped_bin_worker(
    bin_idx: int,
    global_bin_idx: int,
    experiment_uri: str,
    grouped_bin: GroupedBin,
    output_path: Path,
    value_filter: str | None,
    obs_columns: list[str] | None,
    var_columns: list[str] | None,
    var_joinids: list[int] | None,
    x_layer: str,
    output_format: Literal["zarr", "h5ad"],
) -> tuple[int, int, str]:
    """
    Worker function for parallel grouped bin extraction.

    :param bin_idx: Internal index of the bin being processed
    :param global_bin_idx: Global bin index for output file naming
    :param experiment_uri: URI of the SOMA experiment
    :param grouped_bin: GroupedBin to _extract
    :param output_path: Path to write the output file
    :param value_filter: Optional user-provided filter to combine with group filter
    :param obs_columns: Optional list of obs columns to include
    :param var_columns: Optional list of var columns to include
    :param var_joinids: Optional list of var soma_joinids to filter features by
    :param x_layer: Name of the X layer to read
    :param output_format: Output format, either "zarr" or "h5ad"

    :raise SomaExtractError: If extraction fails

    :return: Tuple of (bin_idx, global_bin_idx, output_path_str)
    """
    extract_grouped_bin_to_anndata(
        experiment_uri=experiment_uri,
        grouped_bin=grouped_bin,
        output_path=output_path,
        value_filter=value_filter,
        obs_columns=obs_columns,
        var_columns=var_columns,
        var_joinids=var_joinids,
        x_layer=x_layer,
        output_format=output_format,
    )

    return bin_idx, global_bin_idx, str(output_path)


def extract_grouped_bins(
    *,
    curriculum_metadata: GroupedCurriculumMetadata,
    output_dir: Path,
    partition_index: int = 0,
    max_bins_per_partition: int | None = None,
    output_format: Literal["zarr", "h5ad"] = "h5ad",
    max_workers: int | None = None,
    verbose: bool = True,
) -> None:
    """
    Extract grouped bins in parallel with support for distributed workers.

    Each grouped bin becomes one output file. Cells from the same group stay together.
    No shuffle stage is needed for grouped extraction.

    :param curriculum_metadata: SOMA curriculum metadata with grouped_bins
    :param output_dir: Local directory to save AnnData files
    :param partition_index: Zero-based partition index for distributed workers
    :param max_bins_per_partition: Number of bins per partition (for distribution)
    :param output_format: Output format - "zarr" or "h5ad" (default: "h5ad")
    :param max_workers: Maximum number of parallel workers
    :param verbose: If False, suppress INFO level logging in parallel workers

    :raise SomaExtractError: If SOMA reads fail
    :raise IOError: If file operations fail
    :raise ValueError: If output_format is invalid or grouped_bins is None
    """
    warnings.filterwarnings(action="ignore", category=ImplicitModificationWarning)

    if curriculum_metadata.grouped_bins is None:
        raise ValueError("curriculum_metadata.grouped_bins cannot be None for grouped extraction")

    if max_bins_per_partition is None:
        max_bins_per_partition = len(curriculum_metadata.grouped_bins)

    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Slice bins for this partition
    bin_slice_start, bin_slice_end = utils.get_block_slice(
        total_items=len(curriculum_metadata.grouped_bins),
        partition_index=partition_index,
        block_size=max_bins_per_partition,
    )

    bins_to_process = curriculum_metadata.grouped_bins[bin_slice_start:bin_slice_end]

    if not bins_to_process:
        logger.info(f"No bins to process for partition {partition_index}")
        return

    logger.info(
        f"Extracting {len(bins_to_process)} grouped bins "
        f"(partition {partition_index}, bins [{bin_slice_start}:{bin_slice_end}]) "
        f"using {max_workers} workers"
    )

    # Determine file extension
    file_ext = ".zarr" if output_format == "zarr" else ".h5ad"

    # Use spawn method for multiprocessing
    mp_context = multiprocessing.get_context("spawn")

    completed = 0
    failed = []

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp_context,
        initializer=child_init,
        initargs=(logging.getLevelName(logger.level), verbose),
    ) as executor:
        futures = {}
        for bin_idx, grouped_bin in enumerate(bins_to_process):
            global_bin_idx = bin_slice_start + bin_idx
            output_path = output_dir / f"extract_{global_bin_idx:06d}{file_ext}"

            future = executor.submit(
                _extract_grouped_bin_worker,
                bin_idx,
                global_bin_idx,
                curriculum_metadata.experiment_uri,
                grouped_bin,
                output_path,
                curriculum_metadata.value_filter,
                curriculum_metadata.obs_columns,
                curriculum_metadata.var_columns,
                curriculum_metadata.var_joinids,
                curriculum_metadata.x_layer,
                output_format,
            )
            futures[future] = bin_idx

        # Wait for completion
        for future in tqdm(as_completed(futures), total=len(bins_to_process), desc="Extracting grouped bins"):
            try:
                bin_idx, global_bin_idx, output_path = future.result()
                completed += 1
            except Exception as e:
                bin_idx = futures[future]
                failed.append((bin_idx, str(e)))
                logger.error(f"Grouped bin {bin_idx} failed: {e}")

    if failed:
        error_msg = f"Failed to extract {len(failed)} grouped bins: {failed}"
        logger.error(error_msg)
        raise SomaExtractError(error_msg)

    logger.info(f"Successfully extracted all {len(bins_to_process)} grouped bins to {output_dir}")
