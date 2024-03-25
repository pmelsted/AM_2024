#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
import numpy as np
from os.path import exists
import scanpy as sc
import pandas as pd


from methods import get_top_nn
from methods import get_batch_loc_of_top_nn

from scipy.stats import mannwhitneyu


from scipy.spatial.distance import cdist
import multiprocessing


# %%
def diff_clust(clust, clust_comb):
    """

    Get difference between 2 clusters.

    Parameters
    ----------
    clust : Dataframe column
        Clusters of all cels in order, original.
    clust_comb : Dataframe column
        Clusters of all cels in order, after batch effect.

    Returns
    -------
    TYPE
        Metric that measures how different clustering is.

    """
    clust_arr = clust
    clust_arr = np.array(clust.cat.codes)
    clust_comb_arr = np.array(clust_comb.cat.codes)
    mat = []
    mat_comb = []
    for i in range(clust.shape[0]):
        one_cell = np.zeros(clust_arr.shape[0])
        one_cell_comb = np.zeros(clust_arr.shape[0])
        one_cell[np.where(clust_arr == clust_arr[i])] = 1
        one_cell_comb[np.where(clust_comb_arr == clust_comb_arr[i])] = 1
        mat.append(one_cell)
        mat_comb.append(one_cell_comb)
    mat = np.array(mat)
    mat_comb = np.array(mat_comb)

    cluster_size = []
    cell_cluster_size = []
    for i in range(np.unique(clust).shape[0]):
        cluster_size.append([i, np.where(str(i) == clust)[0].shape[0]])
    cluster_size = np.array(cluster_size)
    for i in range(clust.shape[0]):
        cell_cluster_size.append(cluster_size[clust_arr[i], 1])

    nr_matches = []
    nr_matches_zeros = []
    for i in range(clust.shape[0]):
        ones = np.where(mat[i] == 1)[0]
        zeros = np.where(mat[i] == 0)[0]
        nr_matches.append(np.where(mat_comb[i, ones] == 1)[0].shape[0] / ones.shape[0])
        nr_matches_zeros.append(np.where(mat_comb[i, zeros] == 0)[0].shape[0] / zeros.shape[0])

    return np.array(nr_matches), np.array(nr_matches_zeros), np.array(cell_cluster_size)


def cluster_change(adata, adata_comb):
    """

    Get confusion matrix of cluster differences before and after batch correction.

    Parameters
    ----------
    adata : Anndata
        uncorrected data.
    adata_comb : Anndata
        corrected dataf for a given method.

    Returns
    -------
    p : Dataframe
        Confusion matrix of clusters in long form.


    """

    clust1 = np.array(adata.obs["leiden"].cat.codes)

    clust2 = np.array(adata_comb.obs["leiden"].cat.codes)
    nr_matches3 = pd.DataFrame(list(zip(clust1, clust2)), columns=["clust1", "clust2"])
    nr_matches3["clust1"] = nr_matches3["clust1"].astype("category")
    nr_matches3["clust2"] = nr_matches3["clust2"].astype("category")
    p = nr_matches3.groupby(["clust1", "clust2"]).size()
    return p


def cluster_change_single(adata, adata_single):
    nr_matches, nr_matches_zeros, cell_cluster_size = diff_clust(
        adata[adata_single.obs["orig_cells"], :].obs["leiden"], adata_single.obs["leiden"]
    )
    # adata[adata_single.obs["orig_cells"],:]
    clust1 = np.array(adata[adata_single.obs["orig_cells"], :].obs["leiden"].cat.codes)
    clust2 = np.array(adata_single.obs["leiden"].cat.codes)
    # nr_matches2 = np.column_stack((nr_matches,clust1,clust2))
    nr_matches3 = pd.DataFrame(list(zip(clust1, clust2)), columns=["clust1", "clust2"])
    nr_matches3["clust1"] = nr_matches3["clust1"].astype("category")
    nr_matches3["clust2"] = nr_matches3["clust2"].astype("category")
    p = nr_matches3.groupby(["clust1", "clust2"]).size()
    # nr_matches3.groupby(["clust2"]).count()
    return p, nr_matches, nr_matches_zeros, cell_cluster_size


def get_statistically_diff_genes(adata, clust1, clust2):
    res = []
    for i in np.array(adata.var.index.get_indexer(adata.var.loc[adata.var["highly_variable"] == True].index)):
        test = mannwhitneyu(
            np.array(adata.X.todense()[adata.obs.index.get_indexer(adata.obs.loc[adata.obs["leiden"] == clust1].index), i]),
            np.array(adata.X.todense()[adata.obs.index.get_indexer(adata.obs.loc[adata.obs["leiden"] == clust2].index), i]),
        )

        if not np.isnan(test[0][0]) and test[1][0] * 2000 < 0.05:
            res.append(i)
    return np.array(res)


def get_consensus_clusters(adata_orig, adata_methods):
    """Get the consensus clusters.

    Parameters
    ----------
    adata_orig : AnnData object
        the original adata, unbatched.
    adata_methods : List of AnnData objects
        List of AnnData object, with different batch correction methods.

    Returns
    -------
    arr_methods : List of lists
        Each sublist has 2 numpy arrays. First array is the consensus clustering.
        I.e. the ratio of cells in each cluster that do NOT go into the "consensus cluster".
        Smaller is better
        The second has the size ratio of each cluster of the whole.
        0.1 means that the given cluster has 10% of all observations, i.e. cells

    """
    arr_methods = []
    for j in adata_methods:
        if (
            np.in1d(adata_orig.obs.index.values, j.obs.index.values).nonzero()[0].shape[0] == adata_orig.shape[0]
            and adata_orig.shape[0] == j.shape[0]
        ):
            # adata_norms_2 = cdist(k[:, k.var.highly_variable].X.todense(), j[:, j.var.highly_variable].X.todense())
            cluster_diff = cluster_change(adata_orig, j)
        else:
            if len(adata_orig.obs.index.values[0]) != len(j.obs.index.values[0]):
                match = adata_orig.obs.index.isin(j.obs.index.str[:-2].values)
                match_in_comp = np.in1d(j.obs.index.str[:-2].values, adata_orig.obs.index.values).nonzero()[0]
            else:
                match = adata_orig.obs.index.isin(j.obs.index.values)
                match_in_comp = np.in1d(j.obs.index.values, adata_orig.obs.index.values).nonzero()[0]
            match = np.where(match)[0]
            cluster_diff = cluster_change(adata_orig[match, :], j[match_in_comp, :])
        cl_unstacked = np.array(cluster_diff.unstack())
        # arr_clusters = []
        # arr_clusters = (np.sum(cl_unstacked, axis=1) - np.amax(cl_unstacked, axis=1)) / np.sum(cl_unstacked, axis=1)
        arr_clusters = (np.sum(cl_unstacked, axis=1) - np.amax(cl_unstacked, axis=1)) / np.sum(cl_unstacked, axis=1)
        clust_nmbr = np.sum(cl_unstacked, axis=1) / (np.sum(np.sum(cl_unstacked, axis=1)))

        arr_comb_cluster = (np.sum(cl_unstacked, axis=0) - np.amax(cl_unstacked, axis=0)) / np.sum(cl_unstacked, axis=0)
        comb_clust_nmbr = np.sum(cl_unstacked, axis=0) / (np.sum(np.sum(cl_unstacked, axis=0)))

        mm = (np.average(arr_clusters, weights=clust_nmbr) + np.average(arr_comb_cluster, weights=comb_clust_nmbr)) / 2
        # for i in range(cl_unstacked.shape[0]):
        #     arr_clusters.append((np.sum(cl_unstacked[i,:]) - np.amax(cl_unstacked[i,:]))/np.sum(cl_unstacked[i,:]) )
        arr_methods.append([mm])
    # return np.array(arr_methods)
    return arr_methods


def get_consensus_cluster(adata_orig, adata_method):
    """Get the consensus clusters for a single method.

    Parameters
    ----------
    adata_orig : AnnData object
        the original adata, unbatched.
    adata_methods : AnnData object
        Anndata object.

    Returns
    -------
    arr_methods : Numpy Array
        Consensus clustering
    """
    cluster_diff, nr_matches, nr_matches_zeros, clust_size = cluster_change(adata_orig, adata_method)
    cl_unstacked = np.array(cluster_diff.unstack())
    arr_clusters = []
    for i in range(cl_unstacked.shape[0]):
        arr_clusters.append((np.sum(cl_unstacked[i, :]) - np.amax(cl_unstacked[i, :])) / np.sum(cl_unstacked[i, :]))
    arr_clusters
    return np.array(arr_clusters)


def marker_gene_diff(adata1, adata2, cl1, cl2):
    """.

    Parameters
    ----------
    adata1 : Anndata
        First adata file to compare.
    adata2 : Anndata
        Second adata file.
    cl1 : TYPE
        index of cluster 1.
    cl2 : TYPE
        index of cluster 2.

    Returns
    -------
    genes : list
        list of indices of common genes in 2 adata for a pair of clusters.

    """
    adata1.uns["log1p"]["base"] = None
    adata2.uns["log1p"]["base"] = None
    if "rank_genes_groups" not in adata1.uns:
        sc.tl.rank_genes_groups(adata1, "leiden", method="wilcoxon")
    if "rank_genes_groups" not in adata2.uns:
        sc.tl.rank_genes_groups(adata2, "leiden", method="wilcoxon")
    genes1 = np.array(pd.DataFrame(adata1.uns["rank_genes_groups"]["names"]).head(2000))
    genes2 = np.array(pd.DataFrame(adata2.uns["rank_genes_groups"]["names"]).head(2000))

    genes = []
    for i in range(genes1.shape[0]):
        k = np.where(genes1[i, int(cl1)] == genes2[:, int(cl2)])[0]
        if k.shape[0] != 0:
            genes.append(np.where(genes1[i, int(cl1)] == genes2[:, int(cl2)])[0][0])
        else:
            genes.append(-1)
    return genes


def nn_rank_change(adata, adata_list, neuro=False):
    """.

    Function to get NN rank displacement

    Parameters
    ----------
    adata : Anndata
        Original adata.
    adata_list : TYPE
        List of adata of other methods to compare.
    neuro : boolean, optional
        If working on the mouse brain data or not. The default is False.

    Returns
    -------
    ndarray
        Median rank displacement of all methods.

    """
    mean_rank_change = []
    adata_norms = []  # placeholder
    if neuro:
        if exists("data/neuro_dist.npy"):
            adata_norms = np.load("data/neuro_dist.npy")

            # data.close()
        else:
            adata_norms = cdist(adata[:, adata.var.highly_variable].X.todense(), adata[:, adata.var.highly_variable].X.todense())
            print("calculated adata norms")
            np.save("data/neuro_dist.npy", adata_norms)
    else:
        adata_norms = cdist(adata[:, adata.var.highly_variable].X.todense(), adata[:, adata.var.highly_variable].X.todense())

    top_nn = get_top_nn(adata_norms, 30)
    temp_adata_list = adata_list.copy()

    if multiprocessing.cpu_count() >= 5:
        p = multiprocessing.Pool(5)
    else:
        p = multiprocessing.Pool(1)
    # k = []
    ind = 0
    print("starting pool loop")
    # remove harmony, liger and ligerv2
    if len(adata_list) != 1:
        del temp_adata_list[1]
        del temp_adata_list[2]
        del temp_adata_list[3]
    for result in p.starmap(
        nn_rank_parall_exec,
        zip(temp_adata_list, [adata_norms] * len(temp_adata_list), [adata] * len(temp_adata_list), [top_nn] * len(temp_adata_list)),
    ):
        print(ind)
        mean_rank_change.append(result)
        ind = ind + 1

    # get nn order for bbknn
    # add it first in list
    return np.array(mean_rank_change)


def nn_rank_change_embedding(adata, adata_list):
    """

    Get NN rank change for the embeddings.

    Parameters
    ----------
    adata : Anndata
        Uncorrected adata.
    adata_list : list
        list of other adata files to compare.

    Returns
    -------
    ndarray
        Median NN rank displacement for the embeddings.

    """
    mean_rank_change = []
    adata_emb_norms = cdist(adata.obsm["X_pca"], adata.obsm["X_pca"])
    top_nn = get_top_nn(adata_emb_norms, 30)
    temp_adata_list = adata_list.copy()
    for k in temp_adata_list:
        if np.in1d(adata.obs.index.values, k.obs.index.values).nonzero()[0].shape[0] == adata.shape[0] and adata.shape[0] == k.shape[0]:
            # adata_norms_2 = cdist(k[:, k.var.highly_variable].X.todense(), k[:, k.var.highly_variable].X.todense())
            if "X_harmony" in k.obsm:
                adata_norms_2 = cdist(k.obsm["X_harmony"], k.obsm["X_harmony"])
            elif "X_scVI" in k.obsm:
                adata_norms_2 = cdist(k.obsm["X_scVI"], k.obsm["X_scVI"])
            elif "H_norm" in k.obsm:
                adata_norms_2 = cdist(k.obsm["H_norm"], k.obsm["H_norm"])
            else:
                adata_norms_2 = cdist(k.obsm["X_pca"], k.obsm["X_pca"])
        else:
            if len(adata.obs.index.values[0]) != len(k.obs.index.values[0]):
                match = adata.obs.index.isin(k.obs.index.str[:-2].values)
                match_in_comp = np.in1d(k.obs.index.str[:-2].values, adata.obs.index.values).nonzero()[0]
            else:
                match = adata.obs.index.isin(k.obs.index.values)
                match_in_comp = np.in1d(k.obs.index.values, adata.obs.index.values).nonzero()[0]
            adata_emb_norms = cdist(adata.obsm["X_pca"][match, :], adata.obsm["X_pca"][match, :])

            top_nn = get_top_nn(adata_emb_norms, 30)

            if "X_harmony" in k.obsm:
                adata_norms_2 = cdist(k.obsm["X_harmony"][match_in_comp, :], k.obsm["X_harmony"][match_in_comp, :])
            elif "X_scVI" in k.obsm:
                adata_norms_2 = cdist(k.obsm["X_scVI"][match_in_comp, :], k.obsm["X_scVI"][match_in_comp, :])
            elif "H_norm" in k.obsm:
                adata_norms_2 = cdist(k.obsm["H_norm"][match_in_comp, :], k.obsm["H_norm"][match_in_comp, :])
            else:
                adata_norms_2 = cdist(k.obsm["X_pca"][match_in_comp, :], k.obsm["X_pca"][match_in_comp, :])

        top_nn_new = get_batch_loc_of_top_nn(top_nn, adata_emb_norms, adata_norms_2)
        nn_rank_change = np.zeros(top_nn_new.shape)
        for i in range(top_nn_new.shape[0]):
            for j in range(top_nn_new.shape[1]):
                nn_rank_change[i, j] = np.abs(j - top_nn_new[i, j])
        arr = []
        for i in nn_rank_change:
            arr.append(np.median(i))
        mean_rank_change.append(arr)
    return np.array(mean_rank_change)


def nn_rank_parall_exec(k, adata_norms, adata, top_nn):
    """

    Paralell version of nearest neighbor calculations.

    Parameters
    ----------
    k : Anndata
        original data.
    adata_norms : ndarray
        cell to cell distance matrix.
    adata : Anndata
        adata for single method.
    top_nn : ndarray
        top nearest neighbors.

    Returns
    -------
    arr : ndarray
        Median NN displacement.

    """
    if np.in1d(adata.obs.index.values, k.obs.index.values).nonzero()[0].shape[0] == adata.shape[0] and adata.shape[0] == k.shape[0]:
        adata_norms_2 = cdist(k[:, k.var.highly_variable].X.todense(), k[:, k.var.highly_variable].X.todense())
    else:
        if len(adata.obs.index.values[0]) != len(k.obs.index.values[0]):
            match = adata.obs.index.isin(k.obs.index.str[:-2].values)
            match_in_comp = np.in1d(k.obs.index.str[:-2].values, adata.obs.index.values).nonzero()[0]
        else:
            match = adata.obs.index.isin(k.obs.index.values)
            match_in_comp = np.in1d(k.obs.index.values, adata.obs.index.values).nonzero()[0]
        match = np.where(match)[0]
        adata_norms = cdist(adata[match, adata.var.highly_variable].X.todense(), adata[match, adata.var.highly_variable].X.todense())

        top_nn = get_top_nn(adata_norms, 30)

        adata_norms_2 = cdist(k[match_in_comp, k.var.highly_variable].X.todense(), k[match_in_comp, k.var.highly_variable].X.todense())

    top_nn_new = get_batch_loc_of_top_nn(top_nn, adata_norms, adata_norms_2)
    nn_rank_change = np.zeros(top_nn_new.shape)
    for i in range(top_nn_new.shape[0]):
        for j in range(top_nn_new.shape[1]):
            nn_rank_change[i, j] = np.abs(j - top_nn_new[i, j])
    arr = []
    for i in nn_rank_change:
        arr.append(np.median(i))
    return arr


# write list to binary file
def write_list(ls, filename):
    # store list in binary file so 'wb' mode
    file = []
    if exists(filename):
        # "with" statements are very handy for opening files.
        with open(filename, "rb") as fp:
            file = pickle.load(fp)
    file.append(ls)
    with open(filename, "wb") as fp:
        pickle.dump(file, fp)
        # print('Done writing list into a binary file')


# Read list to memory
def read_list(filename):
    # for reading also binary mode is important
    with open(filename, "rb") as fp:
        n_list = pickle.load(fp)
        return n_list


# %%
def main(inputs, wildcards, outputs):
    """.

    Calculate consensus clusters, NN rank displacement
    and leiden clusters and write to pickle files for each iteration.


    Parameters
    ----------
    inputs : adata files
    wildcards : list
        snakemake wildcards.
    outputs : list
        Snakemake outputs.

    Returns
    -------
    None.

    """
    neuro = False
    # resamp = False
    if "simul-neuro" in wildcards[0]:
        cc_file_string = "data/cc_file_simul_neuro.pickle"
        nn_rank_file_string = "data/nn_rank_file_simul_neuro.pickle"
        nn_rank_emb_file_string = "data/nn_rank_emb_file_simul_neuro.pickle"
        orig_leiden_clust_string = "data/leiden_simul_neuro.csv"
    elif "neuro" in wildcards[0]:
        cc_file_string = "data/cc_file_neuro.pickle"
        nn_rank_file_string = "data/nn_rank_file_neuro.pickle"
        nn_rank_emb_file_string = "data/nn_rank_emb_file_neuro.pickle"
        orig_leiden_clust_string = "data/leiden_neuro.csv"
        neuro = True
        # resamp = True
        # neuro_resamp = True
    elif "simul-pbmc" in wildcards[0]:
        cc_file_string = "data/cc_file_simul_pbmc.pickle"
        nn_rank_file_string = "data/nn_rank_file_simul_pbmc.pickle"
        nn_rank_emb_file_string = "data/nn_rank_emb_file_simul_pbmc.pickle"
        orig_leiden_clust_string = "data/leiden_simul_pbmc.csv"
        # neuro_resamp = True
    elif "pbmc4k" in wildcards[0]:
        cc_file_string = "data/cc_file_pbmc4k.pickle"
        nn_rank_file_string = "data/nn_rank_file_pbmc4k.pickle"
        nn_rank_emb_file_string = "data/nn_rank_emb_file_pbmc4k.pickle"
        orig_leiden_clust_string = "data/leiden_pbmc4k.csv"
    elif "heart" in wildcards[0]:
        cc_file_string = "data/cc_file_heart.pickle"
        nn_rank_file_string = "data/nn_rank_file_heart.pickle"
        nn_rank_emb_file_string = "data/nn_rank_emb_file_heart.pickle"
        orig_leiden_clust_string = "data/leiden_heart.csv"
    else:
        cc_file_string = "data/cc_file.pickle"
        nn_rank_file_string = "data/nn_rank_file.pickle"
        nn_rank_emb_file_string = "data/nn_rank_emb_file.pickle"
        orig_leiden_clust_string = "data/leiden.csv"

    names = []
    files = []
    for i in inputs:
        names.append(i[5:-5])
        files.append(sc.read_h5ad(i))
    # files = [files[0],files[5]]

    # list of files = l
    cons_clust = get_consensus_clusters(files[0], files[1:])
    print("finished consensus clusters")

    write_list(cons_clust, cc_file_string)

    print("saved consensus clusters")

    nn_ranks = nn_rank_change(files[0], files[2:], neuro)
    print("finished ranks calculations - ", wildcards[0])
    nn_ranks_emb = nn_rank_change_embedding(files[0], files[2:])
    print("finished embedding ranks calculations - ", wildcards[0])

    write_list(nn_ranks, nn_rank_file_string)
    write_list(nn_ranks_emb, nn_rank_emb_file_string)
    # write_list(nn_rank_comb, nn_rank_comb_file_string)
    if not exists(orig_leiden_clust_string):
        # pd.DataFrame(data=adata.obs["leiden"]).to_csv(outputs[5])
        files[0].obs["leiden"].to_csv(orig_leiden_clust_string)

    # %%


# %%
main(snakemake.input, snakemake.wildcards, snakemake.output)
with open(snakemake.output[0], "w") as f:
    f.write("Create a new text file!")
