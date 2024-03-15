#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import methods


# %%
def main():
    # print("derp")
    # if(snakemake.params["input"] == "adata"):
    adata = methods.run_pp("heart", outputs=snakemake.output)
    # adata = methods.run(True,outputs = snakemake.output)
    methods.apply_method("none", adata, outputs=snakemake.output)
    adata_single = methods.run_pp_single_batch("heart")
    # adata_batched.write("data/python_to_r.h5ad",as_dense=["X"],force_dense=True)
    # np.savetxt("data/variable_genes.csv", np.array(adata_batched.var.sort_values(by="dispersions", ascending=[False]).index), delimiter=",",fmt='%s')
    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc
        import numpy as np
        import pandas as pd

        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)

    adata.write(snakemake.output[0])
    adata_single.write(snakemake.output[1])


if __name__ == "__main__":
    main()
