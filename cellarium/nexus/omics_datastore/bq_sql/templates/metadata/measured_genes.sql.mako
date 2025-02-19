select  original_feature_id,
        cas_feature_index,
        index
from `${project}.${dataset}.${extract_table_prefix}__extract_feature_info`
order by index 