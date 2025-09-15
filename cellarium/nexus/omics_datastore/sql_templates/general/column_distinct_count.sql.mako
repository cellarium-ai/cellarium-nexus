<%!
    from nexus.omics_datastore.bq_ops import constants
%>
-- Count distinct values of a column in a base table
select count(distinct(${column_name})) as distinct_count
from `${project}.${dataset}.${table_name}`
where ${column_name} is not null
