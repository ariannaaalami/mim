import anndata as ad
import mudata as md
import pandas as pd
import pytest

import mim


def test_package_has_version():
    assert mim.__version__ is not None


def test_jaccard_adata():
    results_df = pd.read_csv("/Users/2021ariannaa/GitHub/mim/tests/jaccard_test_results.csv")
    results_df = results_df.drop("cell_ID.1", axis="columns")
    test_adata = ad.read_h5ad("/Users/2021ariannaa/GitHub/mim/tests/test_data.h5ad")
    jaccard = mim.metrics.jaccard_similarity()
    actual_results_adata = jaccard(test_adata, "cell_ID", "modality", ["ATAC", "GEX"], "dopp")
    actual_results_df = actual_results_adata.obs[["jaccard_similarity"]].drop_duplicates(keep="first")
    actual_results_df["cell_ID"] = actual_results_df.index
    comp_df = actual_results_df.merge(results_df, on="cell_ID", how="inner", suffixes=["_test", "_truth"])
    assert comp_df["jaccard_similarity_test"] == comp_df["jaccard_similarity_truth"]


def test_jaccard_mudata():
    results_df = pd.read_csv("/Users/2021ariannaa/GitHub/mim/tests/jaccard_test_results.csv")
    results_df = results_df.drop("cell_ID.1", axis="columns")
    test_mudata = md.read_h5mu("/Users/2021ariannaa/GitHub/mim/tests/test_mudata.h5mu")
    jaccard = mim.metrics.jaccard_similarity()
    actual_results_mudata = jaccard(test_mudata, "cell_ID", "modality", ["ATAC", "GEX"], "dopp")
    actual_results_df = actual_results_mudata.mod["ATAC"].obs[["jaccard_similarity"]]
    actual_results_dopps = actual_results_mudata.mod["ATAC"].obs[actual_results_mudata.mod["ATAC"].obs["dopp"]]
    actual_results_cell_IDs = actual_results_dopps["cell_IDs"]
    actual_results_df["cell_ID"] = actual_results_cell_IDs
    comp_df = actual_results_df.merge(results_df, on="cell_ID", how="inner", suffixes=["_test", "_truth"])
    assert comp_df["jaccard_similarity_test"] == comp_df["jaccard_similarity_truth"]


@pytest.mark.skip(reason="This decorator should be removed when test passes.")
def test_example():
    assert 1 == 0  # This test is designed to fail.
