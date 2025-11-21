import pytest

from cellarium.nexus.omics_datastore import soma_ops


def test_soma_operation_error_is_exception() -> None:
    """
    Verify SomaOperationError is a base Exception.
    """
    assert issubclass(soma_ops.SomaOperationError, Exception)

    # Can be raised and caught
    with pytest.raises(soma_ops.SomaOperationError, match="test error"):
        raise soma_ops.SomaOperationError("test error")


def test_soma_read_error_inherits_from_operation_error() -> None:
    """
    Verify SomaReadError inherits from SomaOperationError.
    """
    assert issubclass(soma_ops.SomaReadError, soma_ops.SomaOperationError)

    with pytest.raises(soma_ops.SomaReadError, match="read failed"):
        raise soma_ops.SomaReadError("read failed")

    # Can be caught as SomaOperationError
    with pytest.raises(soma_ops.SomaOperationError):
        raise soma_ops.SomaReadError("read failed")


def test_soma_write_error_inherits_from_operation_error() -> None:
    """
    Verify SomaWriteError inherits from SomaOperationError.
    """
    assert issubclass(soma_ops.SomaWriteError, soma_ops.SomaOperationError)

    with pytest.raises(soma_ops.SomaWriteError, match="write failed"):
        raise soma_ops.SomaWriteError("write failed")


def test_soma_filter_error_inherits_from_operation_error() -> None:
    """
    Verify SomaFilterError inherits from SomaOperationError.
    """
    assert issubclass(soma_ops.SomaFilterError, soma_ops.SomaOperationError)

    with pytest.raises(soma_ops.SomaFilterError, match="filter error"):
        raise soma_ops.SomaFilterError("filter error")


def test_soma_extract_error_inherits_from_operation_error() -> None:
    """
    Verify SomaExtractError inherits from SomaOperationError.
    """
    assert issubclass(soma_ops.SomaExtractError, soma_ops.SomaOperationError)

    with pytest.raises(soma_ops.SomaExtractError, match="extract failed"):
        raise soma_ops.SomaExtractError("extract failed")


def test_soma_planning_error_inherits_from_operation_error() -> None:
    """
    Verify SomaPlanningError inherits from SomaOperationError.
    """
    assert issubclass(soma_ops.SomaPlanningError, soma_ops.SomaOperationError)

    with pytest.raises(soma_ops.SomaPlanningError, match="planning failed"):
        raise soma_ops.SomaPlanningError("planning failed")


def test_all_exceptions_exported_in_init() -> None:
    """
    Verify all exception types are exported from soma_ops.__init__.
    """
    assert hasattr(soma_ops, "SomaOperationError")
    assert hasattr(soma_ops, "SomaReadError")
    assert hasattr(soma_ops, "SomaWriteError")
    assert hasattr(soma_ops, "SomaFilterError")
    assert hasattr(soma_ops, "SomaExtractError")
    assert hasattr(soma_ops, "SomaPlanningError")


def test_main_exports_in_init() -> None:
    """
    Verify main classes and functions are exported from soma_ops.__init__.
    """
    assert hasattr(soma_ops, "TileDBSOMADataOperator")
    assert hasattr(soma_ops, "build_soma_value_filter")
    assert hasattr(soma_ops, "plan_soma_extract")
    assert hasattr(soma_ops, "shuffle_extracted_chunks")


def test_all_list_contains_expected_names() -> None:
    """
    Verify __all__ contains expected public API names.
    """
    expected_names = {
        "TileDBSOMADataOperator",
        "build_soma_value_filter",
        "plan_soma_extract",
        "shuffle_extracted_chunks",
        "SomaOperationError",
        "SomaReadError",
        "SomaWriteError",
        "SomaFilterError",
        "SomaExtractError",
        "SomaPlanningError",
    }

    assert hasattr(soma_ops, "__all__")
    assert set(soma_ops.__all__) == expected_names
