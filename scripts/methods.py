import scanpy as sc


import pandas as pd
from scipy import sparse
import numpy as np
from scipy.spatial import distance_matrix

from scipy.stats import mannwhitneyu

from scipy.stats import wilcoxon



#%%


def ratio_nn_multiple(match_arr):
    arr = []
    for sub_arr in match_arr:
        filtered = len(list(filter(lambda x: x >= 0, sub_arr)))
        arr.append(filtered / len(sub_arr))
    return np.array(arr)


def get_neighs_from_norm(data, c, n):
    """

    Parameters
    ----------
    data : ndarray
        distance matrix as output from scanpy.
    c : TYPE
        index of cell.
    n : TYPE
        number of nearest neighbours.

    Returns
    -------
    ndarray
        index of top n nearest neighbors.

    """
    cell_row = data[c, :]
    num_nonzero = len(np.nonzero(cell_row)[0])
    top_neighbors = np.argsort(-np.asarray(cell_row))[:num_nonzero]

    # for i in range(len(cell_row)):
    if n <= num_nonzero:
        arr = top_neighbors[-n:]
    else:
        arr = top_neighbors
    return arr[::-1]
    # np.where(top_neighbors < 6)[0]


def get_neighs_from_similarity(data, c, n):
    cell_row = data[c, :]
    num_nonzero = len(np.nonzero(cell_row)[0])
    top_neighbors = np.argsort(-np.asarray(cell_row)).flatten()[:num_nonzero]

    # for i in range(len(cell_row)):
    if n <= num_nonzero:
        arr = top_neighbors[-n:]
    else:
        arr = top_neighbors
    return arr
    # np.where(top_neighbors < 6)[0]


def match_nn_multiple_norm(dist1, dist2, n1, n2):
    results = []
    for i in range(0, dist1.shape[0]):
        # randint = random.randint(0,adata.X.shape[0]-1)
        nn = get_neighs_from_norm(dist1, i, n1)  # get neighbours for cell i
        top_nn = nn  # [:n] # get top n neighbours of cell i
        # print(top_nn)
        cells = []
        nn_comb = get_neighs_from_norm(dist2, i, n2)
        for j in range(len(top_nn)):  # iterate through each top n_i neighbour
            ind = np.where(nn_comb == top_nn[j])
            if len(ind[0]) >= 1:
                cells.append(ind[0][0])
            else:
                cells.append(-1)
            # cells.append(ind[0])

        # if(len(ind[0]) == 1):
        results.append(cells)
        # else:
        #    results.append(-1)
    return results


# get k nearest neighbors for each cell using distance matrix
def getNN(norm, k):
    n = norm.shape[0]
    nn_arr = np.zeros((n, n))
    for i in range(n):
        for j in np.argsort(norm[i, :])[:k]:
            nn_arr[i, j] = norm[i, j]
    return nn_arr


# get total distance from cell c to n - NN
def get_norm_cell_to_nn(cell, n, data, norms):
    # dist = data
    # arr = []
    # for i in range(len(cell)):
    #     arr.append(dist[i,get_neighs(adata,i,n)])
    dist = getNN(norms, n + 1)
    if type(cell) is int:
        return np.sum(dist[cell, get_neighs_from_norm(data, cell, n)])
    else:
        arr = []
        for i in range(len(cell)):
            arr.append(np.sum(data[cell[i], get_neighs_from_norm(data, cell[i], n)]))
        return np.array(arr)


def get_top_nn_norm(norms):
    nn_dist = getNN(norms, 101)
    l = []
    for i in range(nn_dist.shape[0]):
        l.append(nn_dist[i, get_neighs_from_norm(nn_dist, i, 100)])
    return np.array(l)


def get_top_nn_similarity(norms, nn_dist):
    # nn_dist = getNN(norms,101)
    l = []
    for i in range(nn_dist.shape[0]):
        l.append(norms[i, get_neighs_from_similarity(norms, i, 100)])
    return np.array(l)


def get_top_nn(norms, n=None):
    nn_dist = getNN(norms, 101)
    l = []
    for i in range(nn_dist.shape[0]):
        if n is None:
            l.append(get_neighs_from_norm(nn_dist, i, 100))
        else:
            l.append(get_neighs_from_norm(nn_dist, i, n))
    return np.array(l)


# @njit
def get_batch_loc_of_top_nn(top_nn, norms, norms_comb, adata_comb=None):
    all_cells = []
    for i in range(top_nn.shape[0]):
        # cell 1 top NN
        # get vector with n NN for cell c with distance matrix.
        top = get_neighs_from_norm(norms_comb, i, norms_comb.shape[0])
        one_cell = []
        # print(top)
        # iterate over the n NN
        # find where original NN are in the batched NN
        for j in top_nn[i, :]:
            one_cell.append(np.where(top == j)[0][0])
        all_cells.append(one_cell)

    return np.array(all_cells)


def get_orig_place_nn(norms, norms_comb, adata_comb):
    cells_same_batch = []
    cells_diff_batch = []
    batch0 = np.where(adata_comb.obs["batch"] == 0)[0]

    top_nn = get_top_nn(norms)

    for i in range(top_nn.shape[0]):
        # cell 1 top NN

        current_cell_batch = adata_comb.obs["batch"].iloc[i]
        top_same = []
        top_diff = []

        # create 2 lists, those in the top NN that are in same batch, those in diff batch
        # MIGHT BE DIFFERENT SIZE
        #
        for j in top_nn[i, :]:
            if np.where(batch0 == j)[0].shape[0] > 0:  # cell is in batch0
                if current_cell_batch == 0:
                    top_same.append(j)
                else:
                    top_diff.append(j)
            else:  # cell is in batch1
                if current_cell_batch == 1:
                    top_same.append(j)
                else:
                    top_diff.append(j)

        one_cell_same_batch = []
        one_cell_diff_batch = []

        # find what rank they are in non batched nn
        # unoptimized but simpler
        for j in top_same:
            one_cell_same_batch.append(np.where(j == top_nn[i, :])[0][0])

        for j in top_diff:
            one_cell_diff_batch.append(np.where(j == top_nn[i, :])[0][0])

        # for j in top_nn[i,:]:

        cells_same_batch.append(one_cell_same_batch)
        cells_diff_batch.append(one_cell_diff_batch)

    return np.array(cells_same_batch), np.array(cells_diff_batch)


#
def apply_batch(adata, fix, subsample, batch_nr):
    """

    Create pseudo batches from adata.

    Parameters
    ----------
    adata : Anndata
        adata file with data.
    fix : Boolean
        Set random seed or not.
    subsample : Boolean
        Sumsample or not.
    batch_nr : int
        Nr of batches.

    Returns
    -------
    adata : Anndata
        adata with added pseudo batches.

    """
    split = np.zeros((2))
    if fix:
        np.random.seed(42)
    indices = np.random.permutation(adata.X.shape[0])
    ind = np.vstack((np.arange((adata.X.shape[0])), np.zeros(adata.X.shape[0]))).T
    if batch_nr > 2:
        arrs = np.array_split(indices, batch_nr)

        for i in range(batch_nr):
            ind[arrs[i], 1] = i

    else:
        if subsample is True:
            split[0] = 1 / 2
            split[1] = 1 / 4

        else:
            split[0] = 1 / 2
            split[1] = 1 / 2

        a = int(adata.X.shape[0] * split[0])
        b = int(adata.X.shape[0] * split[1])
        split_1, split_2 = indices[:b], indices[a:]

        ind[split_1, 1] = 0
        ind[split_2, 1] = 1
    # adata_comb = adata.copy()
    adata.obs["batch"] = ind[:, 1]
    return adata


def run_pp(data, subsample=False, outputs=None, fix=False, batch_nr=2, diffexp=False):
    """Preprocess adata for given data.

    Parameters.
    ----------
    data : String
        Which data to load.
    subsample : TYPE, optional
        Input for batch function, should one batch be subsampled. The default is False.
    outputs : List of Strings, optional
        Snakemake output strings. The default is None.
    fix : TYPE, optional
        Input for batch function, fix the random seed. The default is False.
    batch_nr : TYPE, optional
        Input for batch function, nr of batches. The default is 2.

    Returns
    -------
    adata : AnnData
        Preprocessed AnnData object.

    """
    min_genes = 200
    min_cells = 3

    if data == "pbmc_orig":
        adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        mt_string = "MT-"
        n_genes_max_filt = 2500
        pct_mt_filt = 5
    elif data == "neuro":
        adata = sc.read_10x_h5("data/1M_neurons_neuron20k.h5")
        min_genes = 500
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata.var_names_make_unique()
        n_genes_max_filt = 6000
        mt_string = "mt-"
        pct_mt_filt = 10
    elif data == "heart":
        output_str = "data/heart/"
        adata = sc.read(output_str + "/output.mtx", cache=False)  # transpose the data*
        adata.var_names = pd.read_csv(
            output_str + "/output.genes.txt",
            header=None,
        )[0]
        adata.obs_names = pd.read_csv(output_str + "/output.barcodes.txt", header=None)[0]

        t2g = pd.read_csv("data/mus_t2g.txt", delimiter="\t", index_col=1, header=None)

        temp = []
        for i in range(adata.var.index.shape[0]):
            if type(t2g.loc[adata.var.index[i], 2]) == pd.pandas.core.series.Series:
                temp.append(t2g.loc[adata.var.index[i], 2][0])
            else:
                temp.append(t2g.loc[adata.var.index[i], 2])
        adata.var = adata.var.set_index(np.array(temp))

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        n_genes_max_filt = 8000
        mt_string = "mt-"
        pct_mt_filt = 65
    elif data == "pbmc_4k":
        output_str = "data/pbmc4k/"
        adata = sc.read(output_str + "/output.mtx", cache=False)  # transpose the data*
        adata.var_names = pd.read_csv(
            output_str + "/output.genes.txt",
            header=None,
        )[0]
        adata.obs_names = pd.read_csv(output_str + "/output.barcodes.txt", header=None)[0]

        t2g = pd.read_csv("data/homo_t2g.txt", delimiter="\t", index_col=1, header=None)

        temp = []
        for i in range(adata.var.index.shape[0]):
            if type(t2g.loc[adata.var.index[i], 2]) == pd.pandas.core.series.Series:
                temp.append(t2g.loc[adata.var.index[i], 2][0])
            else:
                temp.append(t2g.loc[adata.var.index[i], 2])
        adata.var = adata.var.set_index(np.array(temp))

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        n_genes_max_filt = 4000
        mt_string = "MT-"
        pct_mt_filt = 10
    elif data == "jejunum":
        adata = sc.read_10x_h5("data/5k_human_jejunum_CNIK_3pv3_raw_feature_bc_matrix.h5")
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        if "feature_types" in adata.var:
            adata.var.drop(columns=["feature_types"], inplace=True)
        if "genome" in adata.var:
            adata.var.drop(columns=["genome"], inplace=True)
        mt_string = "MT-"
        n_genes_max_filt = 2500
        pct_mt_filt = 10
        n_genes_max_filt = 5000
    elif data == "simul_neuro":
        adata = sc.read_10x_h5("data/1M_neurons_neuron20k.h5")
        min_genes = 500
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata.var_names_make_unique()
        n_genes_max_filt = 6000
        mt_string = "mt-"
        pct_mt_filt = 10
    elif data == "simul_pbmc":
        adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        mt_string = "MT-"
        n_genes_max_filt = 2500
        pct_mt_filt = 5

    if "downsample" in data:
        pass
    else:
        adata.var["mt"] = adata.var_names.str.startswith(mt_string)  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        sc.pp.filter_cells(adata, max_genes=n_genes_max_filt)
        adata = adata[adata.obs.pct_counts_mt < pct_mt_filt, :]
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
    adata = apply_batch(adata, fix, subsample, batch_nr)

    if data == "simul_pbmc":
        top_genes = np.loadtxt("data/top_genes_pbmc_orig.csv", delimiter=",", dtype=str, skiprows=1)
        top_genes = top_genes[:, 3]
        # adata[adata.obs["batch"] == 1, top_genes[:, 3]].X*0.5
        # for i in range(top_genes.shape[0]):
        X = adata[adata.obs["batch"] == 1, top_genes].X.todense() * 0.5
        adata[adata.obs["batch"] == 1, top_genes].X = sparse.csr_matrix(np.around(X, decimals=0))

        # lower expr of topp 100 expressed genes in b cells and t cells
        # in one batch
    if data == "simul_neuro":
        # = np.loadtxt("data/top_genes_pbmc_orig.csv", delimiter=",", dtype=str, skiprows=1)
        # temp preprocessing to get clusters
        sc.pp.filter_cells(adata, min_genes=1)
        sc.pp.filter_genes(adata, min_cells=1)
        adata_temp = adata.copy()
        sc.pp.normalize_total(adata_temp, target_sum=1e4)
        sc.pp.log1p(adata_temp)
        adata_temp.raw = adata_temp
        adata_temp = adata_temp[:, adata_temp.var.highly_variable]
        sc.pp.regress_out(adata_temp, ["total_counts", "pct_counts_mt"])
        sc.pp.scale(adata_temp, max_value=10)
        sc.tl.pca(adata_temp, use_highly_variable=True)
        sc.pp.neighbors(adata_temp, random_state=42)
        sc.tl.leiden(adata_temp, random_state=42)

        sc.tl.rank_genes_groups(adata_temp, "leiden")
        genes = pd.DataFrame(adata_temp.uns["rank_genes_groups"]["names"]).head(100)

        cl1 = np.where(genes == "Ctss")[1][0]
        cl2 = np.where(genes == "C1qb")[1][0]

        top_genes = np.array(genes.iloc[:, cl1])
        print("cluster is ", cl1)
        if cl1 != cl2:
            print("gene mismatch!")
        # Microglia

        X = adata[adata.obs["batch"] == 1, top_genes].X.todense() * 0.5
        adata[adata.obs["batch"] == 1, top_genes].X = sparse.csr_matrix(np.around(X, decimals=0))

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    pd.DataFrame(data=adata.layers["counts"].todense(), index=adata.obs_names, columns=adata.var_names).to_csv(outputs[3])

    adata = adata[:, adata.var.highly_variable]

    if "neuro" in data:
        sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
        # sc.pp.highly_variable_genes(adata,n_top_genes=5000)
        # sc.pp.highly_variable_genes(adata,flavor="seurat_v3",n_top_genes=2000,layer ="counts")
        if not diffexp:
            sc.pp.scale(adata, max_value=10)

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    # sc.tl.pca(adata, use_highly_variable=True)
    sc.tl.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata, random_state=42)
    sc.tl.leiden(adata, random_state=42)

    return adata


def run_pp_single_batch(data):
    """Preprocesses the data the same as the run_pp func but returns only an single batch as an AnnData file.

    Parameters
    ----------
    data : String
        Which data to load.

    Returns
    -------
    adata_single : AnnData object
        Preprocessed AnnData object

    """
    min_genes = 200
    min_cells = 3
    if data == "pbmc_orig":
        adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        mt_string = "MT-"
        n_genes_max_filt = 2500
        pct_mt_filt = 5
    elif data == "neuro":
        adata = sc.read_10x_h5("data/1M_neurons_neuron20k.h5")
        min_genes = 500
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata.var_names_make_unique()
        n_genes_max_filt = 6000
        mt_string = "mt-"
        pct_mt_filt = 10
    elif data == "heart":
        output_str = "data/pbmc4k/"
        adata = sc.read(output_str + "/output.mtx", cache=False)  # transpose the data*
        adata.var_names = pd.read_csv(
            output_str + "/output.genes.txt",
            header=None,
        )[0]
        adata.obs_names = pd.read_csv(output_str + "/output.barcodes.txt", header=None)[0]
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        n_genes_max_filt = 8000
        mt_string = "mt-"
        pct_mt_filt = 65
    elif data == "pbmc_4k":
        output_str = "data/pbmc4k/"
        adata = sc.read(output_str + "/output.mtx", cache=False)  # transpose the data*
        adata.var_names = pd.read_csv(
            output_str + "/output.genes.txt",
            header=None,
        )[0]
        adata.obs_names = pd.read_csv(output_str + "/output.barcodes.txt", header=None)[0]
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        n_genes_max_filt = 3500
        mt_string = "MT-"
        pct_mt_filt = 10
    elif data == "jejunum":
        adata = sc.read_10x_h5("data/5k_human_jejunum_CNIK_3pv3_raw_feature_bc_matrix.h5")
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        mt_string = "MT-"
        n_genes_max_filt = 2500
        pct_mt_filt = 20
        n_genes_max_filt = 5000
    elif data == "resample_heart":
        adata = sc.read_h5ad("data/adata_resampled_heart.h5ad")
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        n_genes_max_filt = 8000
        mt_string = "mt-"
        pct_mt_filt = 65
    elif data == "resample_pbmc":
        adata = sc.read_h5ad("data/adata_resampled_pbmc.h5ad")
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        n_genes_max_filt = 3500
        mt_string = "MT-"
        pct_mt_filt = 65
    elif data == "downsample_heart":
        adata = sc.read_h5ad("data/adata_down_heart.h5ad")
        sc.pp.filter_genes(adata, min_cells=min_cells)
        sc.pp.filter_cells(adata, min_genes=10)
    elif data == "downsample_pbmc":
        adata = sc.read_h5ad("data/adata_down_pbmc.h5ad")
        sc.pp.filter_genes(adata, min_cells=min_cells)
        sc.pp.filter_cells(adata, min_genes=10)
    elif data == "simul_neuro":
        adata = sc.read_10x_h5("data/1M_neurons_neuron20k.h5")
        min_genes = 500
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        adata.var_names_make_unique()
        n_genes_max_filt = 6000
        mt_string = "mt-"
        pct_mt_filt = 10
    elif data == "simul_pbmc":
        adata = sc.read_10x_mtx("data/filtered_gene_bc_matrices/hg19/", var_names="gene_symbols", cache=False)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        mt_string = "MT-"
        n_genes_max_filt = 2500
        pct_mt_filt = 5

    adata.var["mt"] = adata.var_names.str.startswith(mt_string)  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, max_genes=n_genes_max_filt)
    adata = adata[adata.obs.pct_counts_mt < pct_mt_filt, :]
    # adata.layers["counts"] = adata.X.copy()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    np.random.seed(42)

    indices = np.random.permutation(adata.X.shape[0])

    split_1, split_2 = indices[: int(adata.X.shape[0] / 2)], indices[int(adata.X.shape[0] / 2) :]

    adata_single = adata[split_1, :].copy()

    adata_single.obs["orig_cells"] = split_1

    sc.pp.highly_variable_genes(adata_single, flavor="seurat_v3", n_top_genes=2000)
    sc.pp.normalize_total(adata_single, target_sum=1e4)
    sc.pp.log1p(adata_single)

    adata_single = adata_single[:, adata_single.var.highly_variable]
    if data == "neuro":
        sc.pp.regress_out(adata_single, ["total_counts", "pct_counts_mt"])
        # sc.pp.highly_variable_genes(adata,n_top_genes=5000)
        # sc.pp.highly_variable_genes(adata,flavor="seurat_v3",n_top_genes=2000,layer ="counts")

        sc.pp.scale(adata_single, max_value=10)

    if not sparse.issparse(adata_single.X):
        adata_single.X = sparse.csr_matrix(adata_single.X)
    sc.tl.pca(adata_single, use_highly_variable=True)
    sc.pp.neighbors(adata_single, random_state=42)
    sc.tl.leiden(adata_single, random_state=42)
    return adata_single


def apply_method(meth, adata, inputs=None, outputs=None, diffexp=False):
    """
    Run a given batch correction method on the data.

    Also includes the default none method. Which is the initial uncorrected workflow.

    Parameters
    ----------
    meth : string
        name of method.
    adata : Anndata
        data in anndata form.
    inputs : list, optional
        Input from snakemake workflow. The default is None.
    outputs : list, optional
        Outputs from snakemake workflow. The default is None.

    Returns
    -------
    Anndata
        adata file that has been batch corrected.

    """
    adata_comb = adata.copy()

    if meth == "none":
        #output X 
        pd.DataFrame(data=adata.X.todense(), index=adata.obs_names, columns=adata.var_names).to_csv(outputs[2])
        #output PCA
        pd.DataFrame(data=adata.obsm["X_pca"], index=adata.obs_names).to_csv(outputs[5])
        pd.DataFrame(data=adata.layers["counts"].todense(), index=adata.obs_names, columns=adata.var_names).to_csv(outputs[7])
        np.savetxt(outputs[6], np.array(adata_comb.obs["batch"]), delimiter=",", fmt="%s")
        np.savetxt(
            outputs[4],
            np.array(adata_comb.var.sort_values(by="highly_variable_rank", ascending=[True]).index[:2000]),
            delimiter=",",
            fmt="%s",
        )
        # sc.tl.pca(adata_comb)
        return adata_comb

    if meth == "bbknn":
        import bbknn

        sc.tl.pca(adata_comb)
        bbknn.bbknn(adata_comb)
        sc.tl.leiden(adata_comb, random_state=42)
        return adata_comb

    if meth == "combat":
        t = sc.pp.combat(adata_comb, inplace=False)
        # t[t < 0.00001] = 0
        adata_comb.X = sparse.csr_matrix(t)
    if meth == "combatseq":
        #combatseq wants raw counts
        # adata_comb_temp = read_h5ad("r_to_adata.h5ad")
        adata_comb.X = sparse.csr_matrix(np.array(pd.read_csv(inputs[1], sep=" ", header=None)).astype(float))
        sc.pp.normalize_total(adata_comb, target_sum=1e4)
        sc.pp.log1p(adata_comb)
        # adata_comb.X = sparse.csr_matrix(adata_comb.X)
    if meth == "mnn":
        adata_comb.X = sparse.csr_matrix(np.array(pd.read_csv(inputs[1], sep=" ", header=None)))
        # sc.external.pp.mnn_correct(adata_comb, batch_key="batch")[0][0].X
        #adata_comb.X = sparse.csr_matrix(adata_comb.X)
        # sc.pp.neighbors(adata_comb)
        # sc.tl.leiden(adata_comb)
    if meth == "scvi":
        import scvi

        scvi.settings.seed = 42
        scvi.model.SCVI.setup_anndata(adata_comb, layer="counts", batch_key="batch")
        model = scvi.model.SCVI(adata_comb)

        model.train()

        if diffexp:
            X = model.posterior_predictive_sample(adata_comb, n_samples=1).to_scipy_sparse()
        else:
            X = model.get_normalized_expression()
        latent = model.get_latent_representation()
        adata_comb.obsm["X_scVI"] = latent
        adata_comb.X = sparse.csr_matrix(X)
        sc.pp.neighbors(adata_comb, use_rep="X_scVI", random_state=42)
        sc.tl.leiden(adata_comb, random_state=42)
        return adata_comb

    if meth == "scvi_nobatch":
        import scvi

        scvi.settings.seed = 42
        scvi.model.SCVI.setup_anndata(adata_comb, layer="counts")
        model = scvi.model.SCVI(adata_comb)

        model.train()

        if diffexp:
            X = model.posterior_predictive_sample(adata_comb, n_samples=1)
        else:
            X = model.get_normalized_expression()
        latent = model.get_latent_representation()
        adata_comb.obsm["X_scVI"] = latent
        if not sparse.issparse(X):
            adata_comb.X = sparse.csr_matrix(X)
        sc.pp.neighbors(adata_comb, use_rep="X_scVI", random_state=42)
        sc.tl.leiden(adata_comb, random_state=42)
        return adata_comb

    if meth == "harmony":
        # adata_comb_temp = read_h5ad("r_to_adata.h5ad")
        adata_comb.obsm["X_harmony"] = np.array(pd.read_csv(inputs[1], sep=" ", header=None)).T

        # adata_comb.X = sparse.csr_matrix(adata_comb.X)
        sc.pp.neighbors(adata_comb, use_rep="X_harmony", random_state=42)
        sc.tl.leiden(adata_comb, random_state=42)
        return adata_comb
    if meth == "seuratv4":
        x_seurat = pd.read_csv(inputs[1], sep=" ", header=0)
        genes_out_of_order = np.array(x_seurat.columns.values)

        actual_genes = np.array((adata_comb.var.index))
        temp_orig_genes = []
        temp_new_genes = []

        for i in range(actual_genes.shape[0]):
            temp_orig_genes.append("".join(e for e in actual_genes[i] if e.isalnum()))
        for i in range(genes_out_of_order.shape[0]):
            temp_new_genes.append("".join(e for e in genes_out_of_order[i] if e.isalnum()))

        temp_new_genes = np.array(temp_new_genes)
        temp_orig_genes = np.array(temp_orig_genes)

        # FIND GENES THAT START WITH A NUMBER

        genes_with_numbers = [i for i, si in enumerate(temp_orig_genes) if si.startswith(tuple(str(i) for i in range(10)))]
        for i in genes_with_numbers:
            temp_string = "X" + temp_orig_genes[i]
            if np.where(temp_new_genes == temp_string)[0].shape[0] > 0:
                temp_new_genes[np.where(temp_new_genes == temp_string)] = temp_new_genes[np.where(temp_new_genes == temp_string)][0][1:]

        # np.where(temp_new_genes == "LYZ")[0].shape[0] > 0

        # temp_new_genes[np.where(temp_new_genes == "X7SK2")] = temp_new_genes[np.where(temp_new_genes == "X7SK2")][0][1:]

        list_test = []
        for i in range(genes_out_of_order.shape[0]):
            list_test.append(np.where(temp_new_genes[i] == temp_orig_genes)[0][0])

        list_test = np.array(list_test)

        adata_temp = adata_comb[:, list_test].copy()
        adata_temp.X = sparse.csr_matrix(x_seurat)
        # adata_comb.X = sparse.csr_matrix(np.array(np.array(x_seurat)[:,list_test]))
        sc.tl.pca(adata_temp)
        sc.pp.neighbors(adata_temp)
        sc.tl.leiden(adata_temp, random_state=42)
        return adata_temp

    if meth == "seuratv4_v2":
        x_seurat = pd.read_csv(inputs[1], sep=" ", header=0)
        genes_out_of_order = np.array(x_seurat.columns.values)

        actual_genes = np.array((adata_comb.var.index))

        temp_orig_genes = []
        temp_new_genes = []

        list_test = []

        temp_orig_genes = []
        temp_new_genes = []
        for i in range(actual_genes.shape[0]):
            temp_orig_genes.append("".join(e for e in actual_genes[i] if e.isalnum()))
        for i in range(genes_out_of_order.shape[0]):
            temp_new_genes.append("".join(e for e in genes_out_of_order[i] if e.isalnum()))

        list_test = []
        temp_new_genes = np.array(temp_new_genes)
        temp_orig_genes = np.array(temp_orig_genes)

        genes_with_numbers = [i for i, si in enumerate(temp_orig_genes) if si.startswith(tuple(str(i) for i in range(10)))]
        for i in genes_with_numbers:
            temp_string = "X" + temp_orig_genes[i]
            if np.where(temp_new_genes == temp_string)[0].shape[0] > 0:
                temp_new_genes[np.where(temp_new_genes == temp_string)] = temp_new_genes[np.where(temp_new_genes == temp_string)][0][1:]

        for i in range(genes_out_of_order.shape[0]):
            list_test.append(np.where(temp_new_genes[i] == temp_orig_genes)[0][0])

        adata_temp = adata_comb[:, list_test].copy()
        adata_temp.X = sparse.csr_matrix(x_seurat)
        adata_temp.var["highly_variable"].values[:] = True
        sc.tl.pca(adata_temp)
        sc.pp.neighbors(adata_temp, random_state=42)
        sc.tl.leiden(adata_temp, random_state=42)
        return adata_temp

    if meth == "downsample":
        adata_orig = sc.read_h5ad(inputs[0])
        adata_comb.layers["counts"] = adata_comb.X.copy()
        adata_comb.var_names_make_unique()
        # sc.pp.highly_variable_genes(adata_comb, flavor="seurat_v3", n_top_genes=2000)
        sc.pp.normalize_total(adata_comb, target_sum=1e4)
        sc.pp.log1p(adata_comb)

        sorter = np.argsort(adata_comb.var.index)
        m = sorter[np.searchsorted(adata_comb.var.index, adata_orig.var.index.values, sorter=sorter)]
        adata_comb = adata_comb[:, m].copy()

        adata_comb.raw = adata_comb
    if meth == "resample":
        adata_orig = sc.read_h5ad(inputs[0])

        adata_comb.layers["counts"] = adata_comb.X.copy()

        sc.pp.normalize_total(adata_comb, target_sum=1e4)

        sc.pp.log1p(adata_comb)

        sorter = np.argsort(adata_comb.var.index)
        m = sorter[np.searchsorted(adata_comb.var.index, adata_orig.var.index.values, sorter=sorter)]
        adata_comb = adata_comb[:, m].copy()
        adata_comb.raw = adata_comb
    if meth == "liger":
        import pyliger

        adata_temp = adata_comb.copy()
        adata_temp.X = adata_temp.layers["counts"]
        adata1 = adata_temp[adata_temp.obs["batch"] == 0,].copy()
        adata2 = adata_temp[adata_temp.obs["batch"] == 1,].copy()

        adata1.uns["sample_name"] = "sample1"
        adata1.var.index.name = "genes"
        adata2.var.index.name = "genes"
        adata2.uns["sample_name"] = "sample2"
        adata1.obs.index.name = "cells"
        adata2.obs.index.name = "cells"
        adata_liger = pyliger.create_liger([adata1, adata2], remove_missing=False)

        pyliger.normalize(adata_liger)

        genes = adata_liger.adata_list[1].var.index.to_numpy()

        adata_liger.var_genes = adata_liger.adata_list[1].var.index.to_numpy()
        adata_liger.adata_list[0].uns["var_gene_idx"] = adata_liger.adata_list[1].var.index.isin(genes).nonzero()[0]
        adata_liger.adata_list[1].uns["var_gene_idx"] = adata_liger.adata_list[1].var.index.isin(genes).nonzero()[0]

        # pyliger.select_genes(pbmcs, var_thresh=0.2, do_plot=False)
        pyliger.scale_not_center(adata_liger)
        pyliger.optimize_ALS(adata_liger, k=30, rand_seed=42)
        pyliger.quantile_norm(adata_liger)
        pyliger.leiden_cluster(adata_liger)
        return adata_liger

    if meth == "liger_v2":
        import pyliger

        adata_temp = adata_comb.copy()
        adata_temp.X = adata_temp.layers["counts"]
        adata1 = adata_temp[adata_temp.obs["batch"] == 0,].copy()
        adata2 = adata_temp[adata_temp.obs["batch"] == 1,].copy()

        adata1.uns["sample_name"] = "sample1"
        adata1.var.index.name = "genes"
        adata2.var.index.name = "genes"
        adata2.uns["sample_name"] = "sample2"
        adata1.obs.index.name = "cells"
        adata2.obs.index.name = "cells"
        adata_liger = pyliger.create_liger([adata1, adata2])

        pyliger.normalize(adata_liger)

        pyliger.select_genes(adata_liger, do_plot=False)
        pyliger.scale_not_center(adata_liger)
        pyliger.optimize_ALS(adata_liger, k=30, rand_seed=42)

        pyliger.quantile_norm(adata_liger)
        pyliger.leiden_cluster(adata_liger)
        return adata_liger

    if not sparse.issparse(adata_comb.X):
        adata_comb.X = sparse.csr_matrix(adata_comb.X)

    sc.tl.pca(adata_comb)
    sc.pp.neighbors(adata_comb, random_state=42)
    sc.tl.leiden(adata_comb, random_state=42)
    return adata_comb


def get_mw(norms, adata_comb):
    # nn_dist = getNN(norms,101)

    norms_comb = distance_matrix(adata_comb.X.todense(), adata_comb.X.todense())
    top_nn_norms = get_top_nn_norm(norms)
    top_nn_norms_comb = get_top_nn_norm(norms_comb)
    return mannwhitneyu(top_nn_norms[2, :], top_nn_norms_comb[2, :])


def diff_clust(clust, clust_comb, adata):
    clust = clust.cat.codes
    clust = np.array(clust)
    clust_comb = clust_comb.cat.codes
    clust_comb = np.array(clust_comb)
    mat = []
    mat_comb = []
    for i in range(clust.shape[0]):
        one_cell = np.zeros(clust.shape[0])
        one_cell_comb = np.zeros(clust.shape[0])
        one_cell[np.where(clust == clust[i])] = 1
        one_cell_comb[np.where(clust_comb == clust_comb[i])] = 1
        mat.append(one_cell)
        mat_comb.append(one_cell_comb)
    mat = np.array(mat)
    mat_comb = np.array(mat_comb)

    cluster_size = []
    cell_cluster_size = []
    for i in range(np.unique(clust).shape[0]):
        cluster_size.append([i, np.where(str(i) == adata.obs.leiden)[0].shape[0]])
    cluster_size = np.array(cluster_size)
    for i in range(clust.shape[0]):
        cell_cluster_size.append(cluster_size[clust[i], 1])

    nr_matches = []
    nr_matches_zeros = []
    for i in range(clust.shape[0]):
        ones = np.where(mat[i] == 1)[0]
        zeros = np.where(mat[i] == 0)[0]
        nr_matches.append(np.where(mat_comb[i, ones] == 1)[0].shape[0] / ones.shape[0])
        nr_matches_zeros.append(np.where(mat_comb[i, zeros] == 0)[0].shape[0] / zeros.shape[0])
    # return nr_matches
    return np.array(nr_matches), np.array(nr_matches_zeros), np.array(cell_cluster_size)


def cluster_change(adata, adata_comb):
    nr_matches, nr_matches_zeros, cell_cluster_size = diff_clust(adata.obs["leiden"], adata_comb.obs["leiden"], adata)
    clust1 = np.array(adata.obs["leiden"].cat.codes)
    clust2 = np.array(adata_comb.obs["leiden"].cat.codes)
    # nr_matches2 = np.column_stack((nr_matches,clust1,clust2))
    nr_matches3 = pd.DataFrame(list(zip(clust1, clust2)), columns=["clust1", "clust2"])
    nr_matches3["clust1"] = nr_matches3["clust1"].astype("category")
    nr_matches3["clust2"] = nr_matches3["clust2"].astype("category")
    p = nr_matches3.groupby(["clust1", "clust2"]).size()
    # nr_matches3.groupby(["clust2"]).count()
    return p, nr_matches, nr_matches_zeros, cell_cluster_size
