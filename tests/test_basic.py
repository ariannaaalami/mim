import pytest

import mim


def test_package_has_version():
    assert mim.__version__ is not None


# def jaccard_test():
#     jaccard_results_df = pd.read_csv("/Users/2021ariannaa/GitHub/mim/tests/jaccard_test_results.csv")
#     test_adata = ad.read_h5ad("/Users/2021ariannaa/GitHub/mim/tests/test_data.h5ad")
#     jaccard = mim.metrics.jaccard_similarity()
#     jaccard_test_results = jaccard(test_adata, 'cell_ID', 'modality', ['ATAC', 'GEX'])


@pytest.mark.skip(reason="This decorator should be removed when test passes.")
def test_example():
    assert 1 == 0  # This test is designed to fail.
