"""Tests for SOMA validation and ingest pipeline submission utilities."""

from __future__ import annotations

import pytest

from cellarium.nexus.backend.ingest_management.utils import workflows_utils


@pytest.mark.django_db
def test_batch_uris_creates_correct_batches() -> None:
    """
    Ensure _batch_uris creates correct batches.
    """
    uris = [f"gs://bucket/file{i}.h5ad" for i in range(7)]

    # Test with batch size 2: expect 4 batches (2,2,2,1)
    batches = workflows_utils._batch_uris(uris=uris, batch_size=2)
    assert len(batches) == 4
    assert len(batches[0]) == 2
    assert len(batches[1]) == 2
    assert len(batches[2]) == 2
    assert len(batches[3]) == 1


@pytest.mark.django_db
def test_generate_output_uris_creates_correct_paths() -> None:
    """
    Ensure _generate_output_uris generates correct output paths.
    """
    input_uris = [
        "gs://bucket/path/file1.h5ad",
        "gs://bucket/path/file2.h5ad",
        "gs://bucket/path/file3.h5ad",
    ]
    output_dir = "gs://bucket/output/"

    output_uris = workflows_utils._generate_output_uris(input_uris=input_uris, output_dir=output_dir)

    assert len(output_uris) == 3
    assert output_uris[0] == "gs://bucket/output/file1.h5ad"
    assert output_uris[1] == "gs://bucket/output/file2.h5ad"
    assert output_uris[2] == "gs://bucket/output/file3.h5ad"


@pytest.mark.django_db
def test_generate_output_uris_strips_trailing_slash() -> None:
    """
    Ensure _generate_output_uris correctly handles trailing slashes.
    """
    input_uris = ["gs://bucket/file.h5ad"]
    output_dir = "gs://bucket/output///////"

    output_uris = workflows_utils._generate_output_uris(input_uris=input_uris, output_dir=output_dir)

    assert output_uris[0] == "gs://bucket/output/file.h5ad"
