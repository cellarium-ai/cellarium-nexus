"""
Integration tests for SOMA ingest operations.

Test prepare_ingest_plan and ingest_h5ads_partition to verify
end-to-end ingestion works with real SOMA experiments.
"""

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops._ingest.ingest import (
    ingest_h5ads_partition,
    prepare_ingest_plan,
)
from cellarium.nexus.omics_datastore.soma_ops._ingest.preprocessing import (
    sanitize_first_adata_for_schema,
    sanitize_for_ingest,
)
from cellarium.nexus.shared.schemas.omics_datastore import (
    ExperimentVarSchema,
    IngestPlanMetadata,
    IngestSchema,
    ObsSchemaDescriptor,
)

# Tests for prepare_ingest_plan


def test_prepare_ingest_plan_creates_valid_metadata(
    tmp_path,
) -> None:
    """Verify prepare_ingest_plan creates valid IngestPlanMetadata."""
    n_obs = 5
    all_features = ["ENSG0001", "ENSG0002", "ENSG0003"]

    # Create test adata
    X = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0], [0, 0, 5], [1, 1, 0]], dtype=np.int32))
    obs = pd.DataFrame(data={"cell_type": ["a", "b", "c", "d", "e"]}, index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=all_features)  # No columns, only index
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # Sanitize and save h5ad
    sanitize_for_ingest(adata=adata)
    h5ad_path = tmp_path / "test.h5ad"
    adata.write_h5ad(filename=h5ad_path)

    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var),
        x_validation_type="count_matrix",
    )

    experiment_uri = str(tmp_path / "experiment")

    # Prepare first_adata with schema
    first_adata = adata.copy()
    sanitize_first_adata_for_schema(adata=first_adata, ingest_schema=schema)

    plan = prepare_ingest_plan(
        experiment_uri=experiment_uri,
        h5ad_paths=[str(h5ad_path)],
        measurement_name="RNA",
        ingest_schema=schema,
        ingest_batch_size=1,
        first_adata=first_adata,
    )

    assert isinstance(plan, IngestPlanMetadata)
    assert plan.experiment_uri == experiment_uri
    assert plan.total_files == 1
    assert plan.num_partitions == 1
    assert plan.last_partition_size == 1
    assert plan.ingest_batch_size == 1
    assert plan.measurement_name == "RNA"
    assert len(plan.registration_mapping_pickle) > 0

    # Verify experiment was created
    with tiledbsoma.Experiment.open(uri=experiment_uri, mode="r") as exp:
        assert "RNA" in exp.ms


def test_prepare_ingest_plan_computes_correct_partitions_for_multiple_files(
    tmp_path,
) -> None:
    """Verify prepare_ingest_plan computes correct partitions for multiple files."""
    all_features = ["ENSG0001", "ENSG0002", "ENSG0003"]
    h5ad_paths = []

    # Create 5 h5ad files
    for i in range(5):
        X = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0]], dtype=np.int32))
        obs = pd.DataFrame(data={"cell_type": ["a", "b"]}, index=[f"cell_{i}_{j}" for j in range(2)])
        var = pd.DataFrame(index=all_features)  # No columns, only index
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        sanitize_for_ingest(adata=adata)
        h5ad_path = tmp_path / f"file_{i}.h5ad"
        adata.write_h5ad(filename=h5ad_path)
        h5ad_paths.append(str(h5ad_path))

    schema_var_df = pd.DataFrame(index=all_features)  # No columns, only index
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=schema_var_df),
        x_validation_type="count_matrix",
    )

    experiment_uri = str(tmp_path / "experiment")

    # Prepare first_adata with schema
    first_adata = anndata.read_h5ad(filename=h5ad_paths[0])
    sanitize_first_adata_for_schema(adata=first_adata, ingest_schema=schema)

    plan = prepare_ingest_plan(
        experiment_uri=experiment_uri,
        h5ad_paths=h5ad_paths,
        measurement_name="RNA",
        ingest_schema=schema,
        ingest_batch_size=2,
        first_adata=first_adata,
    )

    assert plan.total_files == 5
    assert plan.num_partitions == 3  # ceil(5 / 2) = 3
    assert plan.last_partition_size == 1  # 5 - 2*2 = 1
    assert plan.ingest_batch_size == 2


# Tests for ingest_h5ads_partition


def test_ingest_h5ads_partition_writes_data_to_soma(
    tmp_path,
) -> None:
    """Verify ingest_h5ads_partition writes data to SOMA experiment."""
    all_features = ["ENSG0001", "ENSG0002", "ENSG0003"]

    # Create test adata with known data
    X = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0], [4, 0, 5]], dtype=np.int32))
    obs = pd.DataFrame(
        data={"cell_type": ["type_a", "type_b", "type_a"]},
        index=["cell_0", "cell_1", "cell_2"],
    )
    var = pd.DataFrame(index=all_features)  # No columns, only index
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=var),
        x_validation_type="count_matrix",
    )

    # Sanitize and save h5ad file
    sanitize_for_ingest(adata=adata)
    h5ad_path = tmp_path / "test.h5ad"
    adata.write_h5ad(filename=h5ad_path)

    # Prepare first_adata with schema (expands var to full schema)
    first_adata = adata.copy()
    sanitize_first_adata_for_schema(adata=first_adata, ingest_schema=schema)

    experiment_uri = str(tmp_path / "experiment")

    # Prepare the ingest plan
    plan = prepare_ingest_plan(
        experiment_uri=experiment_uri,
        h5ad_paths=[str(h5ad_path)],
        measurement_name="RNA",
        ingest_schema=schema,
        ingest_batch_size=1,
        first_adata=first_adata,
    )

    # Ingest the partition
    ingest_h5ads_partition(
        ingest_plan=plan,
        partition_index=0,
        local_h5ad_paths=[str(h5ad_path)],
    )

    # Read back and verify data
    with tiledbsoma.Experiment.open(uri=experiment_uri, mode="r") as exp:
        # Check obs
        obs_df = exp.obs.read().concat().to_pandas()
        assert len(obs_df) == 3
        assert "cell_type" in obs_df.columns

        # Check var
        var_df = exp.ms["RNA"].var.read().concat().to_pandas()
        assert len(var_df) == 3

        # Check X data exists
        x_data = exp.ms["RNA"].X["data"]
        assert x_data is not None


def test_ingest_h5ads_partition_appends_data_across_multiple_partitions(
    tmp_path,
) -> None:
    """Verify ingesting multiple partitions correctly appends data."""
    all_features = ["ENSG0001", "ENSG0002"]
    h5ad_paths = []

    schema_var_df = pd.DataFrame(index=all_features)  # No columns, only index
    schema = IngestSchema(
        obs_columns=[ObsSchemaDescriptor(name="cell_type", dtype="str", nullable=False)],
        var_schema=ExperimentVarSchema.from_dataframe(var_df=schema_var_df),
        x_validation_type="count_matrix",
    )

    # Create 3 h5ad files with 2 cells each - sanitize before saving
    for i in range(3):
        X = sp.csr_matrix(np.array([[1, 2], [3, 4]], dtype=np.int32))
        obs = pd.DataFrame(
            data={"cell_type": [f"type_{i}_a", f"type_{i}_b"]},
            index=[f"cell_{i}_0", f"cell_{i}_1"],
        )
        var = pd.DataFrame(index=all_features)  # No columns, only index
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        sanitize_for_ingest(adata=adata)
        h5ad_path = tmp_path / f"file_{i}.h5ad"
        adata.write_h5ad(filename=h5ad_path)
        h5ad_paths.append(str(h5ad_path))

    # Prepare first_adata with schema
    first_adata = anndata.read_h5ad(filename=h5ad_paths[0])
    sanitize_first_adata_for_schema(adata=first_adata, ingest_schema=schema)

    experiment_uri = str(tmp_path / "experiment")

    # Prepare the ingest plan with batch_size=2
    plan = prepare_ingest_plan(
        experiment_uri=experiment_uri,
        h5ad_paths=h5ad_paths,
        measurement_name="RNA",
        ingest_schema=schema,
        ingest_batch_size=2,
        first_adata=first_adata,
    )

    assert plan.num_partitions == 2

    # Ingest partition 0 (files 0 and 1)
    ingest_h5ads_partition(
        ingest_plan=plan,
        partition_index=0,
        local_h5ad_paths=[h5ad_paths[0], h5ad_paths[1]],
    )

    # Ingest partition 1 (file 2)
    ingest_h5ads_partition(
        ingest_plan=plan,
        partition_index=1,
        local_h5ad_paths=[h5ad_paths[2]],
    )

    # Read back and verify all data was ingested
    with tiledbsoma.Experiment.open(uri=experiment_uri, mode="r") as exp:
        obs_df = exp.obs.read().concat().to_pandas()
        # 3 files * 2 cells each = 6 cells
        assert len(obs_df) == 6
