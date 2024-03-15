import anndata

import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_scvi = methods.apply_method("scvi", adata)
    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_scvi.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_scvi_leiden.txt", adata_scvi.obs["leiden"], fmt="%s")
        pd.DataFrame(data=adata_scvi.X.todense(), index=adata_scvi.obs_names, columns=adata_scvi.var_names).to_csv(
            "data/" + snakemake.wildcards[0] + "_scvi_X.csv"
        )
        sc.tl.rank_genes_groups(adata_scvi, "leiden")
        pd.DataFrame(adata_scvi.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_scvi_deg.csv", index=False
        )

    adata_scvi.write(snakemake.output[0])


if __name__ == "__main__":
    main()
