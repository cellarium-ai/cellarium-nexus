<%
    from cellarium.nexus.omics_datastore.bq_ops import constants
    from cellarium.nexus.omics_datastore.bq_ops.bq_sql import mako_helpers as mh

    select_column_names_processed = mh.add_cell_info_required_columns(select_columns)

    if metadata_columns:
        metadata_columns_clauses = [f"json_extract_scalar(c.metadata_extra, '$.{x}') AS {x}" for x in metadata_columns]
        metadata_columns_clause = ", \n".join(metadata_columns_clauses)
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_RANDOMIZED_TABLE_NAME}`
as
${mh.select(select_column_names_processed)},
% if metadata_columns:
    ${metadata_columns_clause},
% endif
% if extract_bin_keys:
    ${mh.build_concat_expression(columns=extract_bin_keys, alias=constants.EXTRACT_BIN_KEY_NAME)},
% endif
    farm_fingerprint(cast(c.id + ${random_seed_offset} as string)) as farm_finger

from `${project}.${dataset}.${constants.BQ_CELL_INFO_TABLE_NAME}` c
join `${project}.${dataset}.${constants.BQ_INGEST_TABLE_NAME}` i on (i.id = c.ingest_id)
${mh.where(filter_statements)}