<%!
    from cellarium.nexus.omics_datastore.bq_ops import constants
%>
<%
    from cellarium.nexus.omics_datastore.bq_sql import mako_helpers as mh

    select_column_names_processed = mh.add_cell_info_required_columns(select_columns)
    select_column_names_processed = mh.remove_leading_alias(select_column_names_processed)
    select_column_names_processed.append("farm_finger")

    if metadata_columns:
        for metadata_column in metadata_columns:
            select_column_names_processed.append(metadata_column)
%>

create or replace table `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_TABLE_NAME}`
partition by range_bucket(extract_bin, generate_array(0, ${partition_bin_count}, ${partition_size}))
cluster by extract_bin
as
% if extract_bin_keys:
with
    numbered as (
        select
            *,
            row_number() over (
                partition by ${constants.EXTRACT_BIN_KEY_NAME}
                order by farm_finger
            ) as rn
        from `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_RANDOMIZED_TABLE_NAME}` c
    ),
    chunked as (
        select
            *,
            cast(floor((rn - 1) / ${extract_bin_size}) as int) as chunk_index_within_category
        from numbered
    ),
    binned as (
        select
            *,
            (dense_rank() over (
                order by ${constants.EXTRACT_BIN_KEY_NAME}, chunk_index_within_category
            ) - 1) as extract_bin
        from chunked
    )
    ${mh.select(select_column_names_processed)},
    extract_bin
from binned
% else:
    ${mh.select(select_column_names_processed)},
    cast(floor((
        row_number() over (
            order by farm_finger
        ) - 1
    ) / ${extract_bin_size}) as int) as extract_bin
from `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_RANDOMIZED_TABLE_NAME}` c
% endif