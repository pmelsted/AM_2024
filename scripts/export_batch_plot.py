#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import scanpy as sc


# %%


# %%
def cluster_change_exp(adata, adata_comb):
    # nr_matches, nr_matches_zeros, cell_cluster_size = diff_clust(adata.obs["leiden"], adata_comb.obs["leiden"], adata)
    clust1 = np.array(adata.obs["leiden"].cat.codes)
    clust2 = np.array(adata_comb.obs["leiden"].cat.codes)
    # nr_matches2 = np.column_stack((nr_matches,clust1,clust2))
    nr_matches3 = pd.DataFrame(list(zip(clust1, clust2)), columns=["clust1", "clust2"])
    nr_matches3["clust1"] = nr_matches3["clust1"].astype("category")
    nr_matches3["clust2"] = nr_matches3["clust2"].astype("category")
    p = nr_matches3.groupby(["clust1", "clust2"]).size()
    # nr_matches3.groupby(["clust2"]).count()
    return p  # , nr_matches, nr_matches_zeros, cell_cluster_size


# %%
def main(neuro=False):
    ls = ["bbknn", "combat","combatseq", "harmony", "liger","liger_v2", "mnn", "scvi", "seurat", "seurat_v2"]
    if neuro:
        ls = ["bbknn", "combat","combatseq", "harmony", "liger", "mnn", "scvi", "seurat"]
    # adata = sc.read_h5ad("data/adata-0.h5ad")
    adata_ls = []
    orig_str = "neuro" if neuro else "pbmc3k"
    adata_orig = sc.read_h5ad(f"data/{orig_str}-0.h5ad")
    cl = []
    for i in ls:
        if i != "":
            tmpstr = "_" + i
        else:
            tmpstr = i
        if neuro:
            adata_ls.append(sc.read_h5ad("data/neuro-0" + tmpstr + ".h5ad"))
        else:
            adata_ls.append(sc.read_h5ad("data/pbmc3k-0" + tmpstr + ".h5ad"))

    if neuro:
        cl_c_str = "data/cluster_concordance_neuro.txt"
    else:
        cl_c_str = "data/cluster_concordance.txt"
    with open(cl_c_str, "w", encoding="utf-8") as f:
        # f.write('Create a new text file!')
        for i,val in enumerate(adata_ls):
            print(adata_ls[i].obs["leiden"].value_counts().shape)
            d = cluster_change_exp(adata_orig, adata_ls[i])
            cl_unstacked = np.array(d.unstack())
            if(neuro):
                np.savetxt(f"data/neuro_{ls[i]}_clustcon.csv", cl_unstacked, delimiter=",", fmt='%d')
            else:
                np.savetxt(f"data/{ls[i]}_clustcon.csv", cl_unstacked, delimiter=",", fmt='%d')
            # cl.append(cl_unstacked)



# %%
if __name__ == "__main__":
    main()
    main(True)


# %%
# create func for plot where ratio of each batch in each square of confusion matrix
def return_ratio(neuro=False):
    ls = ls = ["bbknn", "combat","combatseq", "harmony", "liger", "mnn", "scvi", "seurat"]
    if neuro:
        ls = ["bbknn", "combat","combatseq", "harmony", "liger", "mnn", "scvi", "seurat"]
    # adata = sc.read_h5ad("data/adata-0.h5ad")
    adata_ls = []
    orig_str = "neuro" if neuro else "pbmc3k"
    adata_orig = sc.read_h5ad(f"data/{orig_str}-0.h5ad")
    for i in ls:
        if i != "":
            tmpstr = "_" + i
        else:
            tmpstr = i
        if neuro:
            adata_ls.append(sc.read_h5ad("data/neuro-0" + tmpstr + ".h5ad"))
        else:
            adata_ls.append(sc.read_h5ad("data/pbmc3k-0" + tmpstr + ".h5ad"))
        #if i == "liger":
        #    adata_ls[-1].obs.index = [x[:-2] for x in np.array(adata_ls[-1].obs.index)]

    if neuro:
        cl_c_str = "data/cluster_ratio_neuro.txt"
    else:
        cl_c_str = "data/cluster_ratio.txt"

    # clust = np.array(cluster_change_exp(adata, adata_comb).reset_index())
    
        # f.write('Create a new text file!')
    for i,val in enumerate(adata_ls):
        # sc.tl.umap(adata_ls[i + 1])
        d = cluster_change_exp(adata_orig, adata_ls[i])
        # clust = np.array(cluster_change_exp(adata_ls[0], adata_ls[i + 1]).reset_index())
        print(ls[i])
        d[:] = cluster_ratio(adata_orig, adata_ls[i])[:, 2]

        cl_unstacked = np.array(d.unstack())
        
        if neuro:
            pd.DataFrame(cl_unstacked).to_csv(f"data/neuro_{ls[i]}_clustrat.csv", index=False, header=False)
        else:
            pd.DataFrame(cl_unstacked).to_csv(f"data/{ls[i]}_clustrat.csv", index=False, header=False)
        
        # cl_unstacked = cl_unstacked.T
            # for i in range(cl_unstacked.shape[0]):
            #     for j in range(cl_unstacked.shape[1]):
            #         f.write(str(cl_unstacked[i, j]))
            #         f.write(",")
            #     # f.write(str(cl_unstacked))
            #     f.write("\n")


# %%


def cluster_ratio(adata, adata_comb):
    # i = np.max(clust[:,0])
    # j = np.max(clust[:,1])
    clust = np.array(cluster_change_exp(adata, adata_comb).reset_index())
    c_list = []
    for i in range(np.max(clust[:, 0]) + 1):
        for j in range(np.max(clust[:, 1]) + 1):
            # print(i,j)
            # print(adata[adata.obs["leiden"].cat.codes==i].shape[0])
            # clust1
            size = adata[np.logical_and(adata.obs["leiden"].cat.codes == i, adata_comb.obs["leiden"].cat.codes == j)].obs["batch"].shape[0]
            # print(size)
            mm = (
                adata[np.logical_and(adata.obs["leiden"].cat.codes == i, adata_comb.obs["leiden"].cat.codes == j)].obs["batch"].sum() / size
            )
            # if(size > 9):
            c_list.append([i, j, mm])
    return np.array(c_list)
    # get index of cells


# %%
cluster_ratio(adata,adata_liger)
return_ratio()
