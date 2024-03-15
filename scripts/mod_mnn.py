import anndata
import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])

    adata_mnn = methods.apply_method("mnn", adata, inputs=snakemake.input)
    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_mnn.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_mnn_leiden.txt", adata_mnn.obs["leiden"], fmt="%s")
        pd.DataFrame(data=adata_mnn.X.todense(), index=adata_mnn.obs_names, columns=adata_mnn.var_names).to_csv(
            "data/" + snakemake.wildcards[0] + "_mnn_X.csv"
        )
        sc.tl.rank_genes_groups(adata_mnn, "leiden")
        pd.DataFrame(adata_mnn.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_mnn_deg.csv", index=False
        )

    adata_mnn.write(snakemake.output[0])


if __name__ == "__main__":
    main()
