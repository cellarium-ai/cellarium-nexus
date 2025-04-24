select count(distinct extract_bin) as total_bins,
       (select count(*) 
        from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info`
        where extract_bin = (select max(extract_bin) 
                           from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info`)
       ) as last_bin_size,
       count(*) as total_cells
from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info` 