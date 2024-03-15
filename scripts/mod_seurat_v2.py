
import anndata

import methods

def main():
    
    adata = anndata.read_h5ad(snakemake.input[0])
    adata_seurat = methods.apply_method("seuratv3_v2",adata,inputs=snakemake.input)
    adata_seurat.write(snakemake.output[0])
    
if __name__ == "__main__":
    main()