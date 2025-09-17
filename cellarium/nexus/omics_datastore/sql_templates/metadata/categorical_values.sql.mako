<%!
    from cellarium.nexus.omics_datastore.bq_ops import constants
%>
select distinct ${column_name}
from `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_TABLE_NAME}`
where ${column_name} is not null
order by ${column_name} 