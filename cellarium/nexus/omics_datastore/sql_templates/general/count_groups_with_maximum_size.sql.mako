<%!
    from cellarium.nexus.omics_datastore.bq_ops import constants
%>
<%
    from cellarium.nexus.omics_datastore.bq_sql import mako_helpers as mh


    select_column_names_processed = mh.remove_leading_alias(group_columns)
%>
-- Calculate the total number of groups. Groups are defined by unique combinations of provided columns in
-- `group_columns` argument. Each group will not have more cells than extract_bin_size.
--
-- Parameters:
--   project: GCP project ID
--   dataset: BigQuery dataset name
--   group_columns: Columns used to group data by for counting
--   extract_bin_size: Extract bin size. Each group will not have more cells than this number
--   filter_statements: A dictionary of filters to apply (column_name: value), processed by mh.where

select count(*) as total_groups from (
    with group_counts as (
      ${mh.select(select_column_names_processed)},
        count(*) as n_cells
      from `${project}.${dataset}.${table_name}`
      ${mh.where(filter_statements)}
      ${mh.group_by(select_column_names_processed)}
    )
      ${mh.select(select_column_names_processed)},
      n_cells,
      bin_idx AS bin_number,
      ceil(n_cells / ${extract_bin_size}) AS num_bins,
      (bin_idx - 1) * ${extract_bin_size} + 1 AS bin_start_row,
      least(bin_idx * ${extract_bin_size}, n_cells) AS bin_end_row,
      least(${extract_bin_size}, n_cells - (bin_idx - 1) * ${extract_bin_size}) AS bin_cell_count
    from group_counts,
    unnest(generate_array(1, cast(ceil(n_cells / ${extract_bin_size}) as int64))) as bin_idx
)