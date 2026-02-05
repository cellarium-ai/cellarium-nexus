"""
End-to-end ingest + extract integration test for SOMA operations.
"""

import anndata
import numpy as np
import tiledbsoma

from cellarium.nexus.omics_datastore.soma_ops.data_ingestor import TileDBSOMAIngestor
from cellarium.nexus.omics_datastore.soma_ops.data_operator import TileDBSOMADataOperator


def test_ingest_then_extract_randomized_e2e(
    tmp_path,
    multi_ingest_anndata_inputs: dict[str, object],
) -> None:
    """
    Validate, sanitize, ingest multiple AnnData files, then extract and verify all data.
    """
    ingestor = TileDBSOMAIngestor()
    ingest_schema = multi_ingest_anndata_inputs["ingest_schema"]
    datasets = multi_ingest_anndata_inputs["datasets"]
    expected_gene_sums = multi_ingest_anndata_inputs["expected_gene_sums"]
    expected_total_cells = multi_ingest_anndata_inputs["total_cells"]

    h5ad_paths: list[str] = []
    for adata, h5ad_path in datasets:
        ingestor.validate_and_sanitize_for_ingest(adata=adata, ingest_schema=ingest_schema)
        adata.write_h5ad(filename=h5ad_path)
        h5ad_paths.append(h5ad_path)

    experiment_uri = str(tmp_path / "soma_experiment")

    first_adata = anndata.read_h5ad(filename=h5ad_paths[0])
    plan = ingestor.prepare_ingest_plan(
        experiment_uri=experiment_uri,
        h5ad_paths=h5ad_paths,
        measurement_name="RNA",
        ingest_schema=ingest_schema,
        ingest_batch_size=2,
        first_adata=first_adata,
    )

    ingestor.ingest_h5ads_partition(
        ingest_plan=plan,
        local_h5ad_paths=h5ad_paths[:2],
    )
    ingestor.ingest_h5ads_partition(
        ingest_plan=plan,
        local_h5ad_paths=h5ad_paths[2:],
    )

    with tiledbsoma.Experiment.open(uri=experiment_uri, mode="r") as exp:
        x_layers = list(exp.ms["RNA"].X.keys())
    assert x_layers
    x_layer = "X" if "X" in x_layers else ("data" if "data" in x_layers else x_layers[0])

    operator = TileDBSOMADataOperator(experiment_uri=experiment_uri)
    curriculum = operator.prepare_curriculum_metadata(
        filters=None,
        range_size=50,
        extract_bin_size=32,
        shuffle_ranges=True,
        obs_columns=["cell_type", "tissue", "donor_id"],
        var_columns=["feature_name"],
        x_layer=x_layer,
    )

    output_dir = tmp_path / "extract_output"
    temp_dir = tmp_path / "extract_temp"
    operator.extract_randomized(
        curriculum_metadata=curriculum,
        output_dir=output_dir,
        temp_dir=temp_dir,
        max_workers_extract=1,
        max_workers_shuffle=1,
        cleanup_temp=True,
        verbose=False,
    )

    output_files = sorted(output_dir.glob("extract_*.h5ad"))
    assert len(output_files) == curriculum.num_bins
    assert len(output_files) > 10

    total_cells = 0
    seen_joinids: set[int] = set()
    actual_gene_sums = {name: 0 for name in expected_gene_sums}

    for output_path in output_files:
        chunk = anndata.read_h5ad(filename=output_path)
        total_cells += chunk.n_obs
        seen_joinids.update(int(v) for v in chunk.obs.index)

        assert {"cell_type", "tissue", "donor_id"}.issubset(chunk.obs.columns)
        assert "feature_name" in chunk.var.columns

        col_sums = np.asarray(chunk.X.sum(axis=0)).ravel()
        for feature_name, total in zip(chunk.var["feature_name"].tolist(), col_sums):
            actual_gene_sums[feature_name] += int(total)

    assert total_cells == expected_total_cells
    assert len(seen_joinids) == expected_total_cells
    assert actual_gene_sums == expected_gene_sums
