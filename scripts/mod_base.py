import methods

# import os

import re


def main():
    adata = methods.run_pp("pbmc_orig", outputs=snakemake.output)
    # dd = re.findall(r"\d+", snakemake.wildcards[0])
    # if len(dd) == 2:
    #     import scanpy as sc

    #     ls = [2, 2, 2, 2, 2, 10, 10, 10, 10, 10, 30, 30, 30, 30, 30, 75, 75, 75, 75, 75, 100, 100, 100, 100, 100]
    #     sc.tl.pca(adata, use_highly_variable=True, n_comps=ls[int(dd[1])])
    #     sc.pp.neighbors(adata, random_state=42)
    #     sc.tl.leiden(adata, random_state=42)

    methods.apply_method("none", adata, outputs=snakemake.output)
    adata_single = methods.run_pp_single_batch("pbmc_orig")

    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc
        import numpy as np
        import pandas as pd

        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)
        # sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
        genes = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200)
        clust1 = np.where(genes == "MS4A1")[1][0]
        clust2 = np.where(genes == "CD8A")[1][0]

        adata.uns["diffexp_clusters"] = [str(clust1), str(clust2)]

    adata.write(snakemake.output[0])
    adata_single.write(snakemake.output[1])


if __name__ == "__main__":
    main()
