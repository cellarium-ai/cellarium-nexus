from typing import Final

# Mapping from SOMA/Numpy dtypes to Pandas nullable dtypes
# This ensures that when we create a Series of NA values or validate columns,
# we use a dtype that supports pd.NA (e.g. Int32 instead of int32).
PANDAS_NULLABLE_DTYPES: Final[dict[str, str]] = {
    "bool": "boolean",
    "int8": "Int8",
    "int16": "Int16",
    "int32": "Int32",
    "int64": "Int64",
    "uint8": "UInt8",
    "uint16": "UInt16",
    "uint32": "UInt32",
    "uint64": "UInt64",
    "float32": "Float32",
    "float64": "Float64",
    "str": "string",
    "string": "string",
}
