<%
    from cellarium.nexus.omics_datastore.bq_ops.bq_sql import mako_helpers as mh
%>
-- Calculate the total count of cells based on provided filters.
--
-- Parameters:
--   project: GCP project ID
--   dataset: BigQuery dataset name
--   table_name: Name of the cell_info table
--   filter_statements: A dictionary of filters to apply (column_name: value), processed by mh.where
select
    count(1) as total_cells
from
    `${project}.${dataset}.${table_name}` c
${mh.where(filter_statements)}
;