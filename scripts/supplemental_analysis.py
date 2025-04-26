#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scanpy as sc
import pandas as pd
import numpy as np
import scipy

#process output from R to create a new Seurat preprocessed h5ad
def parse_seurat_data_to_scanpy(orig,seurat_pp_csv,scanpy_pp_h5ad,batch_csv,scanpy_variable_genes_csv,raw_counts,new_pca,new_counts,new_h5ad):

    #Seurat reorders cells and genes so this is function changes the cell order, and more to be able to compare

    #process original to get gene names working
    adata_orig = sc.read_10x_mtx(
        orig, var_names="gene_symbols",
        cache=True)

    adata_orig.var.index
    adata_orig.var["names"] = adata_orig.var.index

    adata2_x = pd.read_csv(seurat_pp_csv, sep=" ")
    adata1 = sc.read_h5ad(scanpy_pp_h5ad)
    adata1_df = adata1.to_df()
    

    adata2_x_reordered = adata2_x.reindex(adata1_df.index)
    


    adata_new = sc.AnnData(adata2_x_reordered,adata2_x_reordered.index.to_frame(), adata2_x_reordered.columns.to_frame())
    #check if cells same order
    print(np.array_equal(np.array(adata1.obs.index),np.array(adata_new.obs.index)))

    column_names = list(adata_new.obs.columns)
    column_names[0] = 'names'  # Change first column
    adata_new.obs.columns = column_names
    adata_new.obs.columns[0]

    
    column_names = list(adata_new.var.columns)
    column_names[0] = 'gene_ids'  # Change first column

    adata_new.var.columns = column_names
    adata_new.obs["batch"] = np.array(pd.read_csv(batch_csv, sep=",",header=None)[0])

    

    adata_new.var["names"] = adata_orig.var[adata_orig.var["gene_ids"].isin(adata_new.var["gene_ids"])]["names"].values
    adata_new.var = adata_new.var.set_index("names")

    np.savetxt(scanpy_variable_genes_csv,np.array(adata_new.var_names), delimiter=",",fmt="%s")
    #add counts
    adata_raw = pd.read_csv(raw_counts, sep=",",index_col=0)
    adata_new.layers["counts"] = adata_raw.loc[adata_new.obs.index,adata_new.var.index]
    adata_new.var["highly_variable"] = True
    sc.tl.pca(adata_new, use_highly_variable=True)

    adata_new_pd = pd.DataFrame(
        adata_new.X.toarray() if scipy.sparse.issparse(adata_new.X) else adata_new.X,
        index=adata_new.obs_names,
        columns=adata_new.var_names
    )
    adata_new_pca_pd = pd.DataFrame(
        adata_new.obsm["X_pca"].toarray() if scipy.sparse.issparse(adata_new.obsm["X_pca"]) else adata_new.obsm["X_pca"],
        index=adata_new.obs_names,
        columns=[f"PC_{i+1}" for i in range(adata_new.obsm["X_pca"].shape[1])]
    )
    np.array_equal(np.array(adata_new_pca_pd.index),np.array(adata1.obs.index))
    #save correct values to file
    adata_new_pca_pd.to_csv(new_pca, sep=",")
    adata_new_pd.to_csv(new_counts, sep=",")

    #now save h5ad


    adata_new.obs["leiden"]  = adata1.obs["leiden"]

    adata_new.X = scipy.sparse.csr_matrix(adata_new.X)

    adata_new.write(new_h5ad)