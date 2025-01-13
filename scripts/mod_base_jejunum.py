import methods

# import os

import re


def main():
    adata = methods.run_pp("jejunum", outputs=snakemake.output)


    methods.apply_method("none", adata, outputs=snakemake.output)
    adata_single = methods.run_pp_single_batch("jejunum")

    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc
        import numpy as np
        import pandas as pd

        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)
        # sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
        genes = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200)
        clust1 = np.where(genes == "MUC2")[1][0] # goblet cells
        clust2 = np.where(genes == "TRPM5")[1][0] #TUFT cells

        adata.uns["diffexp_clusters"] = [str(clust1), str(clust2)]

    adata.write(snakemake.output[0])
    adata_single.write(snakemake.output[1])


if __name__ == "__main__":
    main()
