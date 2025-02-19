<%!
    from nexus.omics_datastore.bq_ops import constants
%>
with tags as (
  select distinct tag
  from `${project}.${dataset}.${constants.BQ_CELL_INFO_TABLE_NAME}`
),
feature_info as (
  select id, symbol, ensemble_id
  from `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_FEATURE_INFO_TABLE_NAME}`
)
select 
  t.tag,
  array_agg(struct(
    f.id as id,
    f.symbol as symbol,
    f.ensemble_id as ensemble_id
  )) as features
from tags t
cross join feature_info f
group by t.tag
order by t.tag