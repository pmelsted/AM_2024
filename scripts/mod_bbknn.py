#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import anndata
import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])

    adata_bbknn = methods.apply_method("bbknn", adata, inputs=snakemake.input)

    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_bbknn.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_bbknn_leiden.txt", adata_bbknn.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata_bbknn, "leiden")
        pd.DataFrame(adata_bbknn.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_bbknn_deg.csv", index=False
        )
    adata_bbknn.write(snakemake.output[0])


if __name__ == "__main__":
    main()
