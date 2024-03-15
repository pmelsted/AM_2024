#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import anndata
import methods
import numpy as np


def main():
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_pyliger = methods.apply_method("liger_v2", adata)

    merged_adata = adata_pyliger.adata_list[0].concatenate(adata_pyliger.adata_list[1])

    cell_names = [x[:-2] for x in np.array(merged_adata.obs.index)]
    sorter = np.argsort(np.array(cell_names))
    order = sorter[np.searchsorted(np.array(cell_names), np.array(adata.obs.index), sorter=sorter)]

    merged_adata = merged_adata[order]
    merged_adata.obs["leiden"] = merged_adata.obs["cluster"].astype("category")
    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd

        # import numpy as np
        import scanpy as sc

        # merged_adata.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_liger_leiden.txt", merged_adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(merged_adata, "leiden")
        pd.DataFrame(merged_adata.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_liger_deg.csv", index=False
        )
    merged_adata.write(snakemake.output[0])


if __name__ == "__main__":
    main()
