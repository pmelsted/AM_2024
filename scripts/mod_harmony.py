import anndata

import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_harmony = methods.apply_method("harmony", adata, inputs=snakemake.input)
    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_harmony.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_harmony_leiden.txt", adata_harmony.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata_harmony, "leiden")
        pd.DataFrame(adata_harmony.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_harmony_deg.csv", index=False
        )
    adata_harmony.write(snakemake.output[0])


if __name__ == "__main__":
    main()
