<%!
    from cellarium.nexus.omics_datastore.bq_ops import constants
%>
-- Fetch distinct values of a column in a base table
select distinct(${column_name}) as v
from `${project}.${dataset}.${table_name}`
where ${column_name} is not null
order by v
