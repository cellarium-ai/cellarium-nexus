class ComparisonOperators:
    """
    Comparison operators used for filtering data.
    """

    EQUAL = "eq"
    IN = "in"
    NOT_EQUAL = "not_eq"
    NOT_IN = "not_in"
    GREATER_THAN = "gt"
    GREATER_THAN_OR_EQUAL = "gte"
    LESS_THAN = "lt"
    LESS_THAN_OR_EQUAL = "lte"

    CURRENTLY_SUPPORTED = [
        EQUAL,
        IN,
        NOT_EQUAL,
        NOT_IN,
        GREATER_THAN,
        GREATER_THAN_OR_EQUAL,
        LESS_THAN,
        LESS_THAN_OR_EQUAL,
    ]


class TemplateDataDictNames:
    """
    Constant names for keys used in a template data dictionary.
    With a help of this class Mako templates get a more structured and consistent arguments for bq_sql query building.
    """

    PROJECT = "project"
    DATASET = "dataset"
    SELECT = "select_columns"
    FILTERS = "filter_statements"
    REQUIRED_TEMPLATE_KWARGS = [PROJECT, DATASET]


CAS_CELL_INFO_REQUIRED_COLUMN_NAMES = ["c.id", "c.ingest_id"]
