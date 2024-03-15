#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import methods
import os
import pandas as pd
import numpy as np


# %%
def main():
    # run KALLISTO HERE

    idx_file = "data/mus_index.idx"
    t2g_file = "data/mus_t2g.txt"
    output_str = "data/heart/"
    adata_write_str = snakemake.output[0]
    bus_out = output_str + "output_sorted_corrected.bus"

    kallisto_cmd = "kallisto bus -o data/heart/ -i data/mus_index.idx -x 10xv3 -t 20 \
        data/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R1_001.fastq.gz data/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R2_001.fastq.gz \
        data/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R1_001.fastq.gz data/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R2_001.fastq.gz"
    commands0 = [
        "wget -nc https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/mus_musculus.tar.gz",
        "tar --skip-old-files -xf mus_musculus.tar.gz",
        # "cd mus_musculus",
        "cp -n mus_musculus/transcriptome.idx data/mus_index.idx",
        "cp -n mus_musculus/transcripts_to_genes.txt data/mus_t2g.txt"
        # "cd .. && rm -d mus_musculus.tar.gz"
    ]

    commands1 = [
        kallisto_cmd,
        "bustools sort -o data/heart/output.s.bus data/heart/output.bus",
        "bustools correct -o data/heart/output.s.c.bus -w ./data/3M-february-2018.txt data/heart/output.s.bus",
        "bustools sort -o " + bus_out + " data/heart/output.s.c.bus"
        # bustools text -o ./data/output_sorted.txt " + bus_out
    ]
    if not os.path.isfile("data/3M-february-2018.txt"):
        os.system(
            "wget -O data/10x_version2_whitelist.txt https://github.com/BUStools/getting_started/releases/download/getting_started/10xv2_whitelist.txt"
        )

    # data directory to work with and index file and t2g file
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

    adata = methods.run_pp("heart", outputs=snakemake.output)

    # adata = methods.run(True,outputs = snakemake.output)
    methods.apply_method("none", adata, outputs=snakemake.output)
    adata_single = methods.run_pp_single_batch("heart")

    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc

        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)

    adata.write(snakemake.output[0])
    adata_single.write(snakemake.output[1])

    if not os.path.isfile("data/heart_cells"):
        pd.DataFrame({"barcode": adata.obs.index}).to_csv("data/heart_cells.txt", index=False)


if __name__ == "__main__":
    main()
