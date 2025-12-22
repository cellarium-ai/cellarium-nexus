<%!
    from cellarium.nexus.omics_datastore.bq_ops import constants
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_MATRIX_COO_TABLE_NAME}`
partition by range_bucket(extract_bin, generate_array(0, ${partition_bin_count}, ${partition_size}))
cluster by extract_bin
as
select
    b.extract_bin,
    m.cell_id,
    array_agg(struct<feature_id int64, raw_counts int64>
        (ef.id, m.raw_count)) as feature_data
from `${project}.${dataset}.${constants.BQ_RAW_COUNT_MATRIX_COO_TABLE_NAME}` m
join `${project}.${dataset}.${constants.BQ_FEATURE_INFO_TABLE_NAME}` fi 
    on (m.feature_id = fi.id)
join `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_FEATURE_INFO_TABLE_NAME}` ef 
    on (fi.ensemble_id = ef.ensemble_id)
join `${project}.${dataset}.${extract_table_prefix}${constants.BQ_EXTRACT_CELL_INFO_TABLE_NAME}` b 
    on (m.cell_id = b.id)
group by 1,2