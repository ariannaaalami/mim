# Imports
import numpy as np
import pandas as pd
from mudata import MuData
from sklearn.metrics import pairwise_distances


def jaccard_similarity(data, cell_ID, modality, modalities, dopp):
    """Calculate Jaccard Similarity metric of integrated di-modal cell data.

    Parameters
    ----------
    data
        The integrated AnnData object. (AnnData)
    cell_ID
        The name of the obs column containing cell IDs. (string)
    modality
        The name of the obs column containing the modality of the data. (string)
    modalities
        List of modality names, length 2. (list of strings)
    dopp
        The name of the obs column containing bools for whether the cell is a doppelgaenger or not. (string)

    Returns
    -------
    Modified AnnData object with additional jaccard_similarity column in obs. All non-dopp cells have NaN values.
    """
    # check if MuData or AnnData
    if type(data) == MuData:
        # extract embedding and make embedding dataframe
        mod_1 = data.mod[modalities[0]]
        mod_1_dopp = mod_1.obs[mod_1.obs["dopp"]]
        mod_1_embedding = mod_1.obsm["embedding"][mod_1.obs.get_indexer(mod_1_dopp.index)]

        mod_2 = data.mod[modalities[1]]
        mod_2_dopp = mod_2.obs[mod_2.obs["dopp"]]
        mod_2_embedding = mod_2.obsm["embedding"][mod_2.obs.get_indexer(mod_2_dopp.index)]

        mod_1_embedding_df = pd.DataFrame(mod_1_embedding, index=mod_1_dopp[cell_ID])
        mod_2_embedding_df = pd.DataFrame(mod_2_embedding, index=mod_2_dopp[cell_ID])
    else:
        # extract embedding and make embedding dataframe
        embedding = data.obsm["embedding"].copy()

        mod_1 = data.obs[data.obs[modality] == modalities[0] & data.obs[dopp]]
        mod_1_embedding = embedding[data.obs.index.get_indexer(mod_1.index)]
        mod_1_embedding_df = pd.DataFrame(mod_1_embedding, index=mod_1[cell_ID])

        mod_2 = data.obs[data.obs[modality] == modalities[1] & data.obs[dopp]]
        mod_2_embedding = embedding[data.obs.index.get_indexer(mod_2.index)]
        mod_2_embedding_df = pd.DataFrame(mod_2_embedding, index=mod_2[cell_ID])

    # Next, make a distance matrix for each modality
    mod_1_dists = pairwise_distances(mod_1_embedding_df, metric="euclidean")
    mod_2_dists = pairwise_distances(mod_2_embedding_df, metric="euclidean")

    # Next, define inner functions:
    def jaccard(n1, n2):
        """Calculate 'overlap' between two groups of k nearest neighbors.

        Parameters
        ----------
        n1
            List of k-nearest neighbor cell IDs from omics layer 1. (list, bool False if not dopp)
        n2
            List of k-nearest neighbor cell IDs from omics layer 2. (list, bool False if not dopp)

        Returns
        -------
        Float value of Jaccard Similarity, np.nan if cell is not a doppelgaenger.
        """
        if not n1 or not n2:
            return np.nan
        intersection = len(list(set(n1).intersection(n2)))
        union = (len(n1) + len(n2)) - intersection
        return float(intersection) / union

    def k_neighbors(cell, omic, k):
        """Find k nearest cells to parameter cells within the same omics layer.

        Parameters
        ----------
        cell
            The cell barcode and batch (unique identifier) of cell
        omic
            the modality of the cell's information
        k
            the number of nearest neighbors to return

        Returns
        -------
        List of cell IDs of the k nearest cells to the parameter cell. Bool False if not doppelgaenger.
        """
        if not data.obs[dopp].loc[cell]:
            return False
        elif omic == modalities[0]:
            index = mod_1_embedding_df.index.get_loc(cell)
            sorted_dists = np.argsort(mod_1_dists[index])
            neighbors = sorted_dists[1 : k + 1]
            return mod_1_embedding_df.iloc[neighbors].index.values.tolist()
        else:
            index = mod_2_embedding_df.index.get_loc(cell)
            sorted_dists = np.argsort(mod_2_dists[index])
            neighbors = sorted_dists[1 : k + 1]
            return mod_2_embedding_df.iloc[neighbors].index.values.tolist()

    # Next.... make a jaccard similarity column in the data!
    if type(data) == MuData:
        data.mod[modalities[0]].obs["jaccard_similarity"] = (
            data.mod[modalities[0]]
            .obs[cell_ID]
            .apply(lambda x: jaccard(k_neighbors(x, modalities[0], 50), k_neighbors(x, modalities[1], 50)))
        )
        data.mod[modalities[1]].obs["jaccard_similarity"] = (
            data.mod[modalities[1]]
            .obs[cell_ID]
            .apply(lambda x: jaccard(k_neighbors(x, modalities[1], 50), k_neighbors(x, modalities[0], 50)))
        )
        data.update()
    else:
        data.obs["jaccard_similarity"] = data.obs[cell_ID].apply(
            lambda x: jaccard(k_neighbors(x, modalities[0], 50), k_neighbors(x, modalities[1], 50))
        )
    return data
