<%
    from cellarium.nexus.omics_datastore.bq_ops import constants
    from cellarium.nexus.omics_datastore.bq_ops.bq_sql import mako_helpers as mh

    select_column_names_processed = mh.add_cell_info_required_columns(select_columns)
    # Remove column aliases because after creating temporary tables all the columns are
    # part of the single table `c`
    select_column_names_processed = mh.remove_leading_alias(select_column_names_processed)
%>
${mh.select(select_column_names_processed)}
from `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_TABLE_NAME}` c
where extract_bin between ${start_bin} and ${end_bin}
order by farm_finger 