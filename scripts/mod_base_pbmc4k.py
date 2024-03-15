#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import methods
import os
import pandas as pd
import numpy as np


# %%
def main():
    # run KALLISTO HERE

    # data directory to work with and index file and t2g file
    output_str = "data/pbmc4k/"
    idx_file = "data/homo_index.idx"
    t2g_file = "data/homo_t2g.txt"

    adata_write_str = snakemake.output[0]
    # final bus output file
    bus_out = output_str + "output_sorted_corrected.bus"
    # run the kallisto bus to create the initial busfile from fastqs
    kallisto_cmd = (
        "kallisto bus -o "
        + output_str
        + " -i data/homo_index.idx -x 10xv2 -t 20 \
        ./data/fastqs/pbmc4k_S1_L001_R1_001.fastq.gz ./data/fastqs/pbmc4k_S1_L001_R2_001.fastq.gz  \
        ./data/fastqs/pbmc4k_S1_L002_R1_001.fastq.gz  ./data/fastqs/pbmc4k_S1_L002_R2_001.fastq.gz"
    )
    # if the transcriptome files and t2g files do not exist get them
    commands0 = [
        "wget -nc https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz",
        "tar --skip-old-files -xf homo_sapiens.tar.gz",
        # "cd mus_musculus",
        "cp -n homo_sapiens/transcriptome.idx data/homo_index.idx",
        "cp -n homo_sapiens/transcripts_to_genes.txt data/homo_t2g.txt"
        # "cd .. && rm -d mus_musculus.tar.gz"
    ]
    # sort correct and sort again, maybe one to many sorts
    commands1 = [
        kallisto_cmd,
        "bustools sort -o " + output_str + "output.s.bus " + output_str + "output.bus",
        "bustools correct -o " + output_str + "output.s.c.bus -w ./data/10x_version2_whitelist.txt " + output_str + "output.s.bus",
        "bustools sort -o " + bus_out + " " + output_str + "output.s.c.bus"
        # "bustools text -o ./data/output_sorted.txt ./data/output_sorted_corrected.bus"
    ]
    # get version 2 whitelist from 10x if it doesnt exist
    if not os.path.isfile("data/10x_version2_whitelist.txt"):
        os.system(
            "wget -O data/10x_version2_whitelist.txt https://github.com/BUStools/getting_started/releases/download/getting_started/10xv2_whitelist.txt"
        )
    if not (os.path.isfile(idx_file) or os.path.isfile(t2g_file)):
        for c in commands0:
            os.system(c)

    if not os.path.isfile(bus_out):
        print("running kallisto cmd")
        for c in commands1:
            os.system(c)
    # run count on the final busfile
    commands2 = [
        "bustools count -o "
        + output_str
        + " -g "
        + t2g_file
        + " -e "
        + output_str
        + "matrix.ec -t "
        + output_str
        + "transcripts.txt  --genecounts "
        + output_str
        + "output.s.c.bus"
    ]

    # exec(open("scripts/sampling_bus.py").read())
    for c in commands2:
        print("running capture and count")
        os.system(c)

    adata = methods.run_pp("pbmc_4k", outputs=snakemake.output)

    # adata = methods.run(True,outputs = snakemake.output)
    methods.apply_method("none", adata, outputs=snakemake.output)
    adata_single = methods.run_pp_single_batch("pbmc_4k")

    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc

        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)

    adata.write(snakemake.output[0])
    adata_single.write(snakemake.output[1])

    if not os.path.isfile("data/pbmc_cells.txt"):
        pd.DataFrame({"barcode": adata.obs.index}).to_csv("data/pbmc_cells.txt", index=False)


if __name__ == "__main__":
    main()
