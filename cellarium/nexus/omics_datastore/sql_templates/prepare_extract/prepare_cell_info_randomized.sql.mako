<%
    from nexus.omics_datastore.bq_ops import constants
    from nexus.omics_datastore.bq_sql import mako_helpers as mh

    select_column_names_processed = mh.add_cell_info_required_columns(select_columns)
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_RANDOMIZED_TABLE_NAME}`
as
${mh.select(select_column_names_processed)},
    farm_fingerprint(cast(c.id + ${random_seed_offset} as string)) as farm_finger
from `${project}.${dataset}.${constants.BQ_CELL_INFO_TABLE_NAME}` c
join `${project}.${dataset}.${constants.BQ_INGEST_TABLE_NAME}` i on (i.id = c.ingest_id)
${mh.where(filter_statements)}