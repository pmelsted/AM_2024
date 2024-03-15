import anndata

import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_combat = methods.apply_method("combat", adata)
    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_combat.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_combat_leiden.txt", adata_combat.obs["leiden"], fmt="%s")
        pd.DataFrame(data=adata_combat.X.todense(), index=adata_combat.obs_names, columns=adata_combat.var_names).to_csv(
            "data/" + snakemake.wildcards[0] + "_combat_X.csv"
        )
        sc.tl.rank_genes_groups(adata_combat, "leiden")
        pd.DataFrame(adata_combat.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_combat_deg.csv", index=False
        )
    adata_combat.write(snakemake.output[0])


if __name__ == "__main__":
    main()
