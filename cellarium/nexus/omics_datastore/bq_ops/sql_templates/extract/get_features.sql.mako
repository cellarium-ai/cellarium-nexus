<%
    from cellarium.nexus.omics_datastore.bq_ops import constants
%>
select 
    id,
    symbol,
    ensemble_id
from `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_FEATURE_INFO_TABLE_NAME}`
order by id 