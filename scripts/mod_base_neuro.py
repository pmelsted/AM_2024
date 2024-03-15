#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# import scanpy as sc
import methods

# import numpy as np
# import os


# %%
def main():
    # print("derp")
    # if(snakemake.params["input"] == "adata"):
    #     adata = methods.run(outputs = snakemake.output)
    # else:
    #     adata = methods.run(True,outputs = snakemake.output)
    adata = methods.run_pp("neuro", outputs=snakemake.output)
    # adata = methods.run(True,outputs = snakemake.output)
    methods.apply_method("none", adata, outputs=snakemake.output)
    adata_single = methods.run_pp_single_batch("neuro")
    # adata_batched.write("data/python_to_r.h5ad",as_dense=["X"],force_dense=True)
    # np.savetxt("data/variable_genes.csv", np.array(adata_batched.var.sort_values(by="dispersions", ascending=[False]).index), delimiter=",",fmt='%s')
    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc
        import numpy as np
        import pandas as pd

        # np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        # sc.tl.rank_genes_groups(adata, "leiden")
        # pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)
        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)
        # sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
        genes = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200)

        clust1 = np.where(genes == "Myl9")[1][0]
        clust2 = np.where(genes == "Ctss")[1][0]

        adata.uns["diffexp_clusters"] = [str(clust1), str(clust2)]
    adata.write(snakemake.output[0])
    adata_single.write(snakemake.output[1])


if __name__ == "__main__":
    main()
