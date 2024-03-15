import anndata

import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_seurat = methods.apply_method("seuratv3", adata, inputs=snakemake.input)

    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_seurat.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_seurat_leiden.txt", adata_seurat.obs["leiden"], fmt="%s")
        pd.DataFrame(data=adata_seurat.X.todense(), index=adata_seurat.obs_names, columns=adata_seurat.var_names).to_csv(
            "data/" + snakemake.wildcards[0] + "_seurat_X.csv"
        )
        sc.tl.rank_genes_groups(adata_seurat, "leiden")
        pd.DataFrame(adata_seurat.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_seurat_deg.csv", index=False
        )

    adata_seurat.write(snakemake.output[0])


if __name__ == "__main__":
    main()
