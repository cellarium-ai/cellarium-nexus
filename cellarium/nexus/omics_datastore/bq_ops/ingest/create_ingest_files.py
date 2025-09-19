"""
Converts input AnnData files to Avro format suitable for loading the Cell Annotation Service Pilot BigQuery schema
version 1.0.0

# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
# https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
# https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py
"""

import concurrent.futures as concurrency
import gc
import json
import logging
import math
import pathlib
import sys
import time
import traceback
import typing as t
from collections.abc import MutableMapping

import anndata
import h5py
import numpy as np
import pandas as pd
import pydantic
from anndata._core.anndata import AnnData
from anndata._io.h5ad import _clean_uns, _read_raw
from anndata._io.specs import read_elem
from fastavro import parse_schema, writer

from cellarium.nexus.omics_datastore.bq_avro_schemas import converter
from cellarium.nexus.omics_datastore.bq_avro_schemas.cell_management import (
    CellInfoBQAvroSchema,
    FeatureInfoBQAvroSchema,
    IngestInfoBQAvroSchema,
)
from cellarium.nexus.omics_datastore.bq_ops import constants, exceptions

# Default value for maximum batch size of avro files when they are being created by avro writer
FLUSH_BATCH_SIZE_DEFAULT = 10000
# Default value for count matrix multiprocessing batch size
COUNT_MATRIX_MULTIPROCESSING_BATCH_SIZE_DEFAULT = 5000
# Default value for feature id lookup.
ORIGINAL_FEATURE_ID_LOOKUP_DEFAULT = "index"


INGEST_METADATA_LIMIT_SIZE = 2**20


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

ParsedFastAvroSchemaType = list[t.Any] | dict[str, t.Any]


class NumpyJSONEncoder(json.JSONEncoder):
    """Special json encoder for numpy types"""

    def default(self, o: t.Any) -> t.Any:
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        else:
            try:
                return json.JSONEncoder.default(self, o)
            except TypeError:
                # In case if the instance is still not serializable just keep track of its type
                return f"Error while encoding type {o.__class__.__name__} while ingest"


def current_milli_time() -> int:
    """
    Return current time in millisecond precision.
    """
    return round(time.time() * 1000)


def get_adata_x_shape(file_path: str) -> tuple[int, int]:
    # This reads minimal metadata and memory-maps the .X
    adata = anndata.read_h5ad(file_path, backed="r")
    shape = adata.shape  # a tuple (n_obs, n_vars)
    adata.file.close()
    return shape


def optimized_read_anndata(input_file_path: pathlib.Path) -> AnnData:
    """
    Read an AnnData object from an `h5ad` file while handling specific attributes.

    This implementation is based on implementation of
    https://github.com/scverse/anndata/blob/6473f2034aa6e28ebc826ceeab15f413b8d294d8/anndata/_io/h5ad.py#L119

    :param input_file_path: File path on a local computer
    :return: AnnData instance
    """
    f = h5py.File(name=input_file_path, mode="r")

    attributes = ["obs", "var", "uns"]
    d = dict(filename=input_file_path, filemode="r")
    d.update({k: read_elem(f[k]) for k in attributes if k in f})

    d["raw"] = _read_raw(f, attrs={"var", "varm"})

    # Backwards compat to <0.7
    if isinstance(f["obs"], h5py.Dataset):
        _clean_uns(d)

    return AnnData(**d)


def optimized_raw_matrix_read_coo(
    input_file_path: pathlib.Path, row_offset: int, end: int
) -> tuple[list[str], list[str], list[int]]:
    """
    Read a sparse matrix from anndata, extract UUIDs for rows and columns.

    :param input_file_path: Path to the anndata file
    :param row_offset: Start index for slicing
    :param end: End index for slicing

    :raise KeyError: If required ID columns are missing

    :return: A tuple containing (row_id_list, col_id_list, data_array)
    """
    adata = optimized_read_anndata(input_file_path=input_file_path)

    # Convert matrix to COO format (sparse matrix format)
    coord = adata.X[row_offset:end, :].tocoo()

    # Ensure the required ID columns exist
    if constants.OBS_NEXUS_ID not in adata.obs.columns:
        raise KeyError(f"Required column '{constants.OBS_NEXUS_ID}' not found in obs data")
    if constants.VAR_NEXUS_ID not in adata.var.columns:
        raise KeyError(f"Required column '{constants.VAR_NEXUS_ID}' not found in var data")

    # Fetch corresponding IDs from metadata
    row_id_list = adata.obs[constants.OBS_NEXUS_ID].iloc[coord.row + row_offset].tolist()
    col_id_list = adata.var[constants.VAR_NEXUS_ID].iloc[coord.col].tolist()

    # Extract data values
    data_array = coord.data.tolist()

    # Cleanup
    adata.file.close()
    del adata
    gc.collect()

    return row_id_list, col_id_list, data_array


def write_avro_with_generator(
    records_generator: t.Callable[[], t.Generator[dict[str, t.Any], t.Any, None]],
    parsed_schema: ParsedFastAvroSchemaType,
    filename: str,
    progress_batch_size: int = 10000,
    flush_batch_size: int = FLUSH_BATCH_SIZE_DEFAULT,
):
    """
    Avro writing worker function. Calls back to `generator` for each row to be written.
    """
    start = current_milli_time()
    with open(filename, "a+b") as out:
        records = []
        counter = 0
        for record in records_generator():
            counter = counter + 1
            records.append(record)

            if counter % flush_batch_size == 0:
                writer(out, parsed_schema, records)
                records = []

            if counter % progress_batch_size == 0:
                end = current_milli_time()
                logger.info(f"    Processed {counter} rows... in {end - start} ms")
                start = end

        if len(records) > 0:
            end = current_milli_time()
            logger.info(f"    Processed {len(records)} rows... in {end - start} ms")
            writer(out, parsed_schema, records)


def write_avro(parsed_schema: ParsedFastAvroSchemaType, output_file_path: str, records: list[dict[str, t.Any]]):
    with open(output_file_path, "a+b") as out:
        writer(fo=out, schema=parsed_schema, records=records)


def write_pandas_to_avro_info_file(
    df: pd.DataFrame,
    df_metadata_extra: pd.DataFrame,
    pydantic_model_class: t.Type[pydantic.BaseModel],
    output_file_name: str,
) -> None:
    """
    Write pandas dataframe to avro file using pydantic model class as a schema

    :param df: Pandas dataframe with all columns from BQ schema
    :param df_metadata_extra: Pandas dataframe with extra metadata for each cell
    :param pydantic_model_class: Schema used for parsing dataframe
    :param output_file_name: Output filename
    """
    avro_schema = converter.pydantic_to_avro(pydantic_model=pydantic_model_class)
    parsed_schema = parse_schema(avro_schema)

    def pandas_row_generator() -> t.Generator[dict[str, t.Any], None, None]:
        for idx, row in df.iterrows():
            yield {**row.to_dict(), constants.OBS_METADATA_EXTRA: df_metadata_extra.loc[idx].dropna().to_json()}

    write_avro_with_generator(
        records_generator=pandas_row_generator, parsed_schema=parsed_schema, filename=output_file_name
    )


def _apply_column_mapping(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """
    Apply column mapping to a DataFrame, handling index mapping as a new column.

    :param df: Input DataFrame to transform
    :param mapping: Dictionary mapping source column names to target names

    :return: Transformed DataFrame with mapped columns
    """
    df = df.copy()

    # Handle index mapping by creating a new column
    if "index" in mapping:
        df[mapping["index"]] = df.index
        mapping = {k: v for k, v in mapping.items() if k != "index"}

    # Rename remaining columns if any
    if mapping:
        df = df.rename(columns=mapping)

    return df


def _required_non_optional_string_fields(model: t.Type[pydantic.BaseModel]) -> set[str]:
    """
    Determine required, non-optional string fields for a Pydantic model.

    A field is considered required and non-optional if its annotation is :class:`str` (not Optional)
    and it has no default value.

    :param model: Pydantic model class to inspect

    :return: Set of field names that must exist and contain strings
    """
    required: set[str] = set()
    for name, field in model.model_fields.items():
        annotation = field.annotation
        if annotation is str and field.is_required():
            required.add(name)
    return required


def _validate_required_string_columns(df: pd.DataFrame, required_cols: set[str], context: str) -> None:
    """
    Validate that required string columns exist and contain non-null string values.

    :param df: Dataframe to validate
    :param required_cols: Set of required string column names
    :param context: Human-readable context for error messages, e.g. "obs" or "var"

    :raise: DataIngestError
    :return: None
    """
    missing = sorted([c for c in required_cols if c not in df.columns])
    if missing:
        raise exceptions.DataValidationError(f"Missing required {context} string columns: {', '.join(missing)}")

    # Check non-null and type
    bad_nulls: list[str] = []
    bad_types: list[str] = []
    for col in required_cols:
        series = df[col]
        if series.isna().any():
            bad_nulls.append(col)
        # Ensure all non-null are strings
        if not series.dropna().map(lambda v: isinstance(v, str)).all():
            bad_types.append(col)

    messages: list[str] = []
    if bad_nulls:
        messages.append(f"columns with null values: {', '.join(sorted(bad_nulls))}")
    if bad_types:
        messages.append(f"columns with non-string values: {', '.join(sorted(bad_types))}")

    if messages:
        raise exceptions.DataValidationError(
            f"Invalid values for required {context} string columns; " + "; ".join(messages)
        )


def _ensure_schema_fields(df: pd.DataFrame, schema_field_names: list[str]) -> pd.DataFrame:
    """
    Ensure all schema fields exist in the dataframe, filling missing columns with None values.

    :param df: Input DataFrame to check and update
    :param schema_field_names: List of required field names from schema

    :return: DataFrame with all required schema fields
    """
    missing_fields = set(schema_field_names) - set(df.columns)
    if missing_fields:
        logger.info(f"Adding missing schema fields with None values: {', '.join(missing_fields)}")
        for field in missing_fields:
            df[field] = pd.Series([pd.NA] * len(df), dtype=pd.Int8Dtype())
    return df


def _process_cell_info_obs(
    adata: AnnData,
    tag: str | None,
    ingest_id: int,
    start_index: int,
    end_index: int,
    column_mapping: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process cell info observations and apply column mapping if provided.

    :param adata: Input anndata object to process
    :param tag: Tag to apply
    :param ingest_id: Ingest ID to apply
    :param start_index: Starting index for cell IDs
    :param end_index: Ending index for cell IDs
    :param column_mapping: Optional mapping of input column names to schema names

    :return: Tuple of (schema_data_df, metadata_extra_df)
    """
    df = adata.obs

    if column_mapping:
        logger.info("Mapping obs columns: " + ", ".join(f"{k} -> {v}" for k, v in column_mapping.items()))
        df = _apply_column_mapping(df=df, mapping=column_mapping)

    schema_field_names = [x for x in CellInfoBQAvroSchema.model_fields if x not in {constants.OBS_METADATA_EXTRA}]
    metadata_extra_columns = list(set(df.columns) - set(schema_field_names))

    # Add Nexus ID, ingest ID, and tag
    df[constants.OBS_NEXUS_ID] = range(start_index, end_index + 1)
    df[constants.OBS_INGEST_ID] = ingest_id
    df[constants.OBS_TAG] = tag

    # Calculate total_mrna_umis as the sum of all counts for each cell
    count_matrix = adata.X
    count_matrix = count_matrix[:]
    total_mrna_umis = count_matrix.sum(axis=1).A1

    df[constants.OBS_TOTAL_MRNA_UMIS] = total_mrna_umis.astype(int)

    logger.info("Ensuring schema fields in obs DataFrame")
    _ensure_schema_fields(df=df, schema_field_names=schema_field_names)

    df_for_schema = df[schema_field_names]
    df_metadata_extra = df[metadata_extra_columns]

    return df_for_schema, df_metadata_extra


def _process_feature_info_var(
    df: pd.DataFrame,
    tag: str | None,
    ingest_id: int,
    start_index: int,
    end_index: int,
    column_mapping: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process feature info variables and apply column mapping if provided.

    :param df: Input DataFrame to process
    :param tag: Tag to apply
    :param ingest_id: Ingest ID to apply
    :param start_index: Starting index for feature IDs
    :param end_index: Ending index for feature IDs
    :param column_mapping: Optional mapping of input column names to schema names

    :return: Tuple of (schema_data_df, metadata_extra_df)
    """
    if column_mapping:
        logger.info("Mapping var columns: " + ", ".join(f"{k} -> {v}" for k, v in column_mapping.items()))
        df = _apply_column_mapping(df=df, mapping=column_mapping)

    schema_field_names = [x for x in FeatureInfoBQAvroSchema.model_fields if x not in {constants.VAR_METADATA_EXTRA}]
    metadata_extra_columns = list(set(df.columns) - set(schema_field_names))

    # Add Nexus ID, ingest ID, and tag
    df[constants.VAR_NEXUS_ID] = range(start_index, end_index + 1)
    df[constants.VAR_INGEST_ID] = ingest_id
    df[constants.VAR_TAG] = tag

    logger.info("Ensuring schema fields in var DataFrame")
    _ensure_schema_fields(df=df, schema_field_names=schema_field_names)

    df_metadata_extra = df[metadata_extra_columns]
    df_for_schema = df[schema_field_names]

    return df_for_schema, df_metadata_extra


def prepare_input_anndata(
    adata: anndata.AnnData,
    tag: str | None,
    ingest_id: int,
    save_path: pathlib.Path,
    cell_index_start: int,
    cell_index_end: int,
    feature_index_start: int,
    feature_index_end: int,
    column_mapping: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Process and prepare AnnData object for ingest.

    :param adata: Input AnnData object
    :param tag: Tag to apply
    :param ingest_id: Ingest ID to apply
    :param save_path: Path to save processed AnnData
    :param cell_index_start: Starting index for cell IDs
    :param cell_index_end: Ending index for cell IDs
    :param feature_index_start: Starting index for feature IDs
    :param feature_index_end: Ending index for feature IDs
    :param column_mapping: Optional mapping of input column names to schema names

    :return: Tuple of processed DataFrames (obs_schema, obs_metadata, var_schema, var_metadata)
    """
    obs_mapping = column_mapping.get("obs_mapping") if column_mapping else None
    var_mapping = column_mapping.get("var_mapping") if column_mapping else None

    df_obs_schema_data, df_obs_metadata_extra = _process_cell_info_obs(
        adata=adata,
        tag=tag,
        ingest_id=ingest_id,
        start_index=cell_index_start,
        end_index=cell_index_end,
        column_mapping=obs_mapping,
    )

    df_var_schema_data, df_var_metadata_extra = _process_feature_info_var(
        df=adata.var,
        tag=tag,
        ingest_id=ingest_id,
        start_index=feature_index_start,
        end_index=feature_index_end,
        column_mapping=var_mapping,
    )

    # Validate required string columns prior to writing AnnData
    obs_required_str = _required_non_optional_string_fields(model=CellInfoBQAvroSchema)
    var_required_str = _required_non_optional_string_fields(model=FeatureInfoBQAvroSchema)

    _validate_required_string_columns(df=df_obs_schema_data, required_cols=obs_required_str, context="obs")
    _validate_required_string_columns(df=df_var_schema_data, required_cols=var_required_str, context="var")

    # Update the AnnData object with processed data
    adata.obs = df_obs_schema_data
    adata.var = df_var_schema_data

    # Save the updated AnnData
    logger.info(f"Saving updated AnnData to {save_path}...")
    adata.write_h5ad(filename=save_path)

    return df_obs_schema_data, df_obs_metadata_extra, df_var_schema_data, df_var_metadata_extra


def dump_cell_info(
    df_obs_schema_data: pd.DataFrame, df_obs_metadata_extra: pd.DataFrame, output_dir: pathlib.Path
) -> None:
    logger.info("Dumping Cell Info into AVRO file...")
    output_file_name = f"{output_dir}/{constants.INGEST_CELL_INFO_FILE_NAME}"
    write_pandas_to_avro_info_file(
        df=df_obs_schema_data,
        df_metadata_extra=df_obs_metadata_extra,
        pydantic_model_class=CellInfoBQAvroSchema,
        output_file_name=output_file_name,
    )


def dump_feature_info(
    df_var_schema_data: pd.DataFrame, df_var_metadata_extra: pd.DataFrame, output_dir: pathlib.Path
) -> None:
    logger.info("Dumping Feature Info into AVRO file...")
    output_file_name = f"{output_dir}/{constants.INGEST_FEATURE_INFO_FILE_NAME}"
    write_pandas_to_avro_info_file(
        df=df_var_schema_data,
        df_metadata_extra=df_var_metadata_extra,
        pydantic_model_class=FeatureInfoBQAvroSchema,
        output_file_name=output_file_name,
    )


def _remove_large_values_from_uns(adata_uns: MutableMapping, uns_key_limit_size_bytes: int) -> None:
    keys_to_replace = []
    for uns_key, uns_value in list(adata_uns.items()):
        value_size = sys.getsizeof(uns_value)
        if value_size > uns_key_limit_size_bytes:
            keys_to_replace.append(uns_key)

    if keys_to_replace:
        logger.info(f"Removing uns keys with large values: {keys_to_replace}")

    for key in keys_to_replace:
        old_value = adata_uns[key]
        value_type = type(old_value)
        value_size = sys.getsizeof(old_value)
        adata_uns[key] = (
            f"Value was removed due to large size. Value type was: {value_type}, value size was: {value_size}."
        )


def _remove_uns_keys(adata_uns: MutableMapping, keys_to_keep: list[str]) -> None:
    keys_to_remove = [key for key in adata_uns.keys() if key not in keys_to_keep]
    logger.info(f"Removing uns keys: {keys_to_remove}")
    for key in keys_to_remove:
        del adata_uns[key]


def dump_ingest_info(
    ingest_id: int,
    adata_uns: MutableMapping,
    output_dir: pathlib.Path,
    metadata_limit_size: int = INGEST_METADATA_LIMIT_SIZE,
    uns_keys_to_keep: list[str] | None = None,
):
    logger.info("Dumping ingest info...")
    ingest_info_obj_to_dump = IngestInfoBQAvroSchema(id=ingest_id)

    output_file_path = f"{output_dir}/{constants.INGEST_INGEST_FILE_NAME}"

    avro_schema = converter.pydantic_to_avro(pydantic_model=IngestInfoBQAvroSchema)
    parsed_schema = parse_schema(avro_schema)

    _remove_large_values_from_uns(adata_uns=adata_uns, uns_key_limit_size_bytes=metadata_limit_size)

    if uns_keys_to_keep is not None:
        _remove_uns_keys(adata_uns=adata_uns, keys_to_keep=uns_keys_to_keep)

    uns_json_dump = json.dumps(obj=adata_uns, cls=NumpyJSONEncoder)

    ingest_info_obj_to_dump.metadata_extra = uns_json_dump

    ingest_info_dict = ingest_info_obj_to_dump.model_dump()
    records = [ingest_info_dict]
    write_avro(parsed_schema=parsed_schema, output_file_path=output_file_path, records=records)


def dump_core_matrix_batch(
    input_file_path: pathlib.Path, row_offset: int, end: int, batch_output_file_name: str, batch_num: int
) -> None:
    """
    Read data, and write the raw X matrix with existing UUID4 identifiers using buffered writing.
    """
    logger.info(f"Starting batch {batch_num} processing rows {row_offset} to {end}")
    row_ids, col_ids, raw_counts = optimized_raw_matrix_read_coo(
        input_file_path=input_file_path, row_offset=row_offset, end=end
    )
    logger.info(f"Batch {batch_num}: Found {len(set(row_ids))} unique cell IDs")
    logger.info(f"Batch {batch_num}: First few cell IDs: {list(set(row_ids))[:5]}")

    start = current_milli_time()
    data_array_int = [int(x) for x in raw_counts]
    logger.info(f"Batch {batch_num} - Convert types... in {current_milli_time() - start} ms")

    # start = current_milli_time()
    buffer = []
    buffer_size = 10000

    with open(batch_output_file_name, "w") as out:
        for row_id, col_id, raw_count in zip(row_ids, col_ids, data_array_int):
            if raw_count != 0:
                buffer.append(f"{row_id},{col_id},{raw_count}\n")

            if len(buffer) >= buffer_size:
                out.writelines(buffer)
                buffer.clear()

        if buffer:
            out.writelines(buffer)

    logger.info(f"Processed Batch {batch_num} in {current_milli_time() - start} ms")


def child_init(log_level: str):
    """
    Called once in each child process at startup.
    Here you configure logging in that worker.
    """
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] [%(processName)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    # logging.info("Child process logging is configured.")


def dump_core_matrix_in_parallel(
    input_file_path: pathlib.Path,
    total_cells: int,
    output_dir: pathlib.Path,
    count_matrix_multiprocessing_batch_size: int = COUNT_MATRIX_MULTIPROCESSING_BATCH_SIZE_DEFAULT,
):
    logger.info("Dumping core matrix into COO format to batched CSV files...")
    num_batches = math.ceil(total_cells / count_matrix_multiprocessing_batch_size)

    # ranges are start-inclusive and end-exclusive
    batches = [
        (
            x,
            x * count_matrix_multiprocessing_batch_size,
            min(total_cells, x * count_matrix_multiprocessing_batch_size + count_matrix_multiprocessing_batch_size),
        )
        for x in range(num_batches)
    ]

    logger.info(f"Total cells: {total_cells}")
    logger.info(f"Batch size: {count_matrix_multiprocessing_batch_size}")
    for batch_num, start, end in batches:
        logger.info(f"Batch {batch_num}: rows {start} to {end}")

    start = current_milli_time()
    logger.info(f"{num_batches} batches will be created to ingest raw count matrix.")

    with concurrency.ProcessPoolExecutor(initializer=child_init, initargs=(logging.INFO,)) as executor:
        futures = []

        for batch_num, i_batch_row_offset, i_batch_end in batches:
            batch_file_name = constants.INGEST_RAW_COUNTS_FILE_NAME_FORMAT.format(batch_index=batch_num)
            batch_output_file_name = f"{output_dir}/{batch_file_name}"
            dump_core_matrix_batch_task_kwargs = {
                "input_file_path": input_file_path,
                "row_offset": i_batch_row_offset,
                "end": i_batch_end,
                "batch_output_file_name": batch_output_file_name,
                "batch_num": batch_num,
            }
            future = executor.submit(dump_core_matrix_batch, **dump_core_matrix_batch_task_kwargs)
            futures.append(future)

        done, not_done = concurrency.wait(futures, return_when=concurrency.ALL_COMPLETED)

        error_raised = False
        for future in done:
            try:
                # Attempt to get the result of the future
                _ = future.result()
            except Exception as e:
                # If an exception is raised, print the exception details
                logger.info(f"Future: {future}")
                logger.info(f"Exception type: {type(e).__name__}")
                logger.info(f"Exception message: {e}")
                # Format and print the full traceback
                traceback.print_exception(type(e), e, e.__traceback__)
                error_raised = True

        if error_raised:
            raise exceptions.DataIngestError("Error occurred during parallel processing of raw count matrix.")

    logger.info(f"    Processed {total_cells} cells... in {current_milli_time() - start} ms")


def create_ingest_files(
    adata_file_path: pathlib.Path,
    tag: str | None,
    cell_info_start_index: int,
    cell_info_end_index: int,
    feature_info_start_index: int,
    feature_info_end_index: int,
    ingest_id: int,
    output_dir: pathlib.Path,
    column_mapping: dict | None = None,
    metadata_limit_size: int = INGEST_METADATA_LIMIT_SIZE,
    uns_keys_to_keep: list[str] | None = None,
) -> dict:
    """
    Create ingest files locally.

    :param adata_file_path: Path to the local .h5ad file
    :param tag: String tag to apply to cell and feature infos
    :param cell_info_start_index: Starting index for cell IDs
    :param cell_info_end_index: Ending index for cell IDs
    :param feature_info_start_index: Starting index for feature IDs
    :param feature_info_end_index: Ending index for feature IDs
    :param ingest_id: The numeric or unique ingest ID
    :param output_dir: Where to store the avro and csv outputs
    :param column_mapping: Optional mapping of input column names to schema names
    :param metadata_limit_size: Max size for uns metadata, etc.
    :param uns_keys_to_keep: List of uns keys to keep in the uns metadata. If None, all keys are kept.

    :raise ValueError: If required columns are missing after mapping

    :return: A dict with pipeline results
    """
    adata = optimized_read_anndata(input_file_path=adata_file_path)

    n_obs = adata.n_obs
    # 2. Data transformations (obs, var, uns)
    logger.info("Adding required columns...")
    adata_updated_path = adata_file_path.parent / "adata-updated.h5ad"
    df_obs_schema_data, df_obs_metadata_extra, df_var_schema_data, df_var_metadata_extra = prepare_input_anndata(
        adata=adata,
        tag=tag,
        ingest_id=ingest_id,
        save_path=adata_updated_path,
        cell_index_start=cell_info_start_index,
        cell_index_end=cell_info_end_index,
        feature_index_start=feature_info_start_index,
        feature_index_end=feature_info_end_index,
        column_mapping=column_mapping,
    )

    # 3. Dump ingest_info, obs, var, matrix
    adata_uns = adata.uns
    adata.file.close()
    del adata

    dump_ingest_info(
        ingest_id=ingest_id,
        adata_uns=adata_uns,
        output_dir=output_dir,
        metadata_limit_size=metadata_limit_size,
        uns_keys_to_keep=uns_keys_to_keep,
    )
    dump_cell_info(
        df_obs_schema_data=df_obs_schema_data, df_obs_metadata_extra=df_obs_metadata_extra, output_dir=output_dir
    )

    dump_feature_info(
        df_var_schema_data=df_var_schema_data, df_var_metadata_extra=df_var_metadata_extra, output_dir=output_dir
    )

    dump_core_matrix_in_parallel(input_file_path=adata_updated_path, total_cells=n_obs, output_dir=output_dir)

    return {
        "output_dir": str(output_dir),
        "files_created": [
            "ingest_info.avro",
            "cell_info.avro",
            "feature_info.avro",
            "raw_counts_*.csv",
        ],
        "adata_uns_clean": json.loads(json.dumps(obj=adata_uns, cls=NumpyJSONEncoder)),
    }
