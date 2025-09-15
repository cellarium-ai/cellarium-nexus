<%!
    from nexus.omics_datastore.bq_ops import constants
%>
-- Count distinct values for multiple columns in a base table
-- Expects: project, dataset, table_name, column_names (list[str])
select
% for i, col in enumerate(column_names):
    count(distinct(${col})) as distinct_${col}${',' if i < len(column_names)-1 else ''}
% endfor
from `${project}.${dataset}.${table_name}`
