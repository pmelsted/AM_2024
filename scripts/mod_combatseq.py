import anndata

import methods


def main():
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_combatseq = methods.apply_method("combatseq", adata, inputs=snakemake.input)
    if "diffexp" in snakemake.wildcards[0]:
        import pandas as pd
        import numpy as np
        import scanpy as sc

        adata_combatseq.uns["log1p"]["base"] = None
        np.savetxt("data/" + snakemake.wildcards[0] + "_combatseq_leiden.txt", adata_combatseq.obs["leiden"], fmt="%s")
        pd.DataFrame(data=adata_combatseq.X.todense(), index=adata_combatseq.obs_names, columns=adata_combatseq.var_names).to_csv(
            "data/" + snakemake.wildcards[0] + "_combatseq_X.csv"
        )
        sc.tl.rank_genes_groups(adata_combatseq, "leiden")
        pd.DataFrame(adata_combatseq.uns["rank_genes_groups"]["names"]).head(200).to_csv(
            "data/" + snakemake.wildcards[0] + "_combatseq_deg.csv", index=False
        )
    adata_combatseq.write(snakemake.output[0])


if __name__ == "__main__":
    main()
