"""
Custom exceptions for SOMA operations.

This module defines SOMA-specific exceptions for better error handling.
"""


class SomaOperationError(Exception):
    """
    Base exception for SOMA operation errors.

    Indicates that a SOMA operation failed for any reason.
    """

    pass


class SomaReadError(SomaOperationError):
    """
    SOMA read operation failed.

    Indicates that reading from SOMA (obs, var, X) encountered an error.
    """

    pass


class SomaWriteError(SomaOperationError):
    """
    SOMA write operation failed.

    Indicates that writing to SOMA encountered an error.
    """

    pass


class SomaFilterError(SomaOperationError):
    """
    Filter translation or application failed.

    Indicates that converting Nexus filters to SOMA value_filter expressions
    or applying filters to SOMA data encountered an error.
    """

    pass


class SomaExtractError(SomaOperationError):
    """
    SOMA data extraction failed.

    Indicates that extracting data from SOMA to AnnData encountered an error.
    """

    pass


class SomaPrepareCurriculumMetadataError(SomaOperationError):
    """
    SOMA extract planning failed.

    Indicates that computing extract plans or joinid ranges encountered an error.
    """

    pass
