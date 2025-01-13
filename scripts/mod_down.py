#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import methods

# import numpy as np
# import os


# %%

import os
import numpy as np
import pandas as pd
import random

# import os
import scanpy as sc

# %%


def sort_cells(adata_new, adata_orig):
    sorter = np.argsort(adata_orig.obs.index)
    m = sorter[np.searchsorted(adata_orig.obs.index, adata_new.obs.index.str[:-2].values, sorter=sorter)]

    va = np.stack((m, np.arange(m.shape[0])), axis=1)
    va = va[va[:, 0].argsort()]
    va[:, 1]
    adata_new = adata_new[va[:, 1]].copy()
    adata_new.obs["batch"] = adata_new.obs.index.str[-1].values.astype("int")
    return adata_new


# %%


def run_sampling(ratio, param):
    # ls = []
    # res = []
    if "heart" in param:
        sampling_depth = 0.5
        output_str = "data/down_heart/"
        idx_file = "data/mus_index.idx"
        t2g_file = "data/mus_t2g.txt"
        bus_out = output_str + "output_sorted_corrected.bus"
        prerun_cells = "data/heart_cells.txt"
        kallisto_cmd = (
            "kallisto bus -o "
            + output_str
            + " -i data/mus_index.idx -x 10xv3 -t 20 \
            data/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R1_001.fastq.gz data/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R2_001.fastq.gz \
            data/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R1_001.fastq.gz data/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R2_001.fastq.gz"
        )
        commands0 = [
            "wget -nc https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/mus_musculus.tar.gz",
            "tar --skip-old-files -xf mus_musculus.tar.gz",
            # "cd mus_musculus",
            # "cp -n mus_musculus/transcriptome.idx data/index.idx",
            # "cp -n mus_musculus/transcripts_to_genes.txt data/t2g.txt"
            "cp -n mus_musculus/transcriptome.idx data/mus_index.idx",
            "cp -n mus_musculus/transcripts_to_genes.txt data/mus_t2g.txt"
            # "cd .. && rm -d mus_musculus.tar.gz"
        ]
        commands1 = [
            kallisto_cmd,
            "bustools sort -o " + output_str + "output.s.bus " + output_str + "output.bus",
            "bustools correct -o " + output_str + "output.s.c.bus -w data/3M-february-2018.txt " + output_str + "output.s.bus",
            "bustools sort -o " + bus_out + " " + output_str + "output.s.c.bus"
            # "bustools text -o ./data/output_sorted.txt ./data/output_sorted_corrected.bus"
        ]

    elif "pbmc4k" in param:
        sampling_depth = 0.5
        output_str = "data/down_pbmc/"
        idx_file = "data/homo_index.idx"
        t2g_file = "data/homo_t2g.txt"

        # adata_write_str = "data/adata_down_pbmc.h5ad"
        bus_out = output_str + "output_sorted_corrected.bus"
        prerun_cells = "data/pbmc_cells.txt"
        kallisto_cmd = (
            "kallisto bus -o "
            + output_str
            + " -i data/homo_index.idx -x 10xv2 -t 20 \
            ./data/fastqs/pbmc4k_S1_L001_R1_001.fastq.gz ./data/fastqs/pbmc4k_S1_L001_R2_001.fastq.gz  \
            ./data/fastqs/pbmc4k_S1_L002_R1_001.fastq.gz  ./data/fastqs/pbmc4k_S1_L002_R2_001.fastq.gz"
        )
        commands0 = [
            "wget -nc https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz",
            "tar --skip-old-files -xf homo_sapiens.tar.gz",
            # "cd mus_musculus",
            "cp -n homo_sapiens/transcriptome.idx data/homo_index.idx",
            "cp -n homo_sapiens/transcripts_to_genes.txt data/homo_t2g.txt"
            # "cd .. && rm -d mus_musculus.tar.gz"
        ]
        commands1 = [
            kallisto_cmd,
            "bustools sort -o " + output_str + "output.s.bus " + output_str + "output.bus",
            "bustools correct -o " + output_str + "output.s.c.bus -w ./data/10x_version2_whitelist.txt " + output_str + "output.s.bus",
            "bustools sort -o " + bus_out + " " + output_str + "output.s.c.bus"
            # "bustools text -o ./data/output_sorted.txt ./data/output_sorted_corrected.bus"
        ]

    if not os.path.isfile("data/10x_version2_whitelist.txt"):
        os.system(
            "wget -O data/10x_version2_whitelist.txt https://github.com/BUStools/getting_started/releases/download/getting_started/10xv2_whitelist.txt"
        )
    if (os.path.isfile(idx_file)) and (os.path.isfile(t2g_file)):
        pass
    else:
        for c in commands0:
            os.system(c)

    if not os.path.isfile(bus_out):
        print("running kallisto cmd")
        for c in commands1:
            os.system(c)

    data = pd.read_csv(prerun_cells)
    data = data["barcode"]
    # batch_index = np.sort(np.random.choice(data.shape[0], size=int(0.5 * data.shape[0]), replace=False))
    batch = np.loadtxt("data/" + param + "_batch.csv").astype(int)

    batch_index = np.where(batch == 0)

    df_a = data.loc[batch_index]
    mask = np.ones(data.shape[0], dtype=bool)
    mask[batch_index] = False
    df_b = data[mask]
    pd.DataFrame(df_a).to_csv(output_str + "cells_a.txt", index=False, header=False)
    pd.DataFrame(df_b).to_csv(output_str + "cells_b.txt", index=False, header=False)

    commands2 = [
        "bustools capture -o " + output_str + "output_a.bus -c " + output_str + "cells_a.txt -b " + bus_out,
        "bustools capture -o " + output_str + "output_b.bus -c " + output_str + "cells_b.txt -b " + bus_out,
        "bustools count -o "
        + output_str
        + "output_b/ -g "
        + t2g_file
        + " -e "
        + output_str
        + "matrix.ec -t "
        + output_str
        + "transcripts.txt  --genecounts "
        + output_str
        + "output_b.bus",
        "bustools count -o "
        + output_str
        + "output_a/ -g "
        + t2g_file
        + " -e "
        + output_str
        + "matrix.ec -t "
        + output_str
        + "transcripts.txt  --genecounts --downsample "
        + str(sampling_depth)
        + " "
        + output_str
        + "output_a.bus",
    ]

    # exec(open("scripts/sampling_bus.py").read())
    for c in commands2:
        print("running capture and count")
        os.system(c)

    adata_a = sc.read(output_str + "output_a/output.mtx", cache=False)  # transpose the dataa*
    adata_a.var_names = pd.read_csv(
        output_str + "output_a/output.genes.txt",
        header=None,
    )[0]
    adata_a.obs_names = pd.read_csv(output_str + "output_a/output.barcodes.txt", header=None)[0]

    adata_b = sc.read(output_str + "output_b/output.mtx", cache=False)  # transpose the data
    adata_b.var_names = pd.read_csv(
        output_str + "output_b/output.genes.txt",
        header=None,
    )[0]
    adata_b.obs_names = pd.read_csv(output_str + "output_b/output.barcodes.txt", header=None)[0]

    adata = adata_a.concatenate(adata_b)
    adata_orig = sc.read_h5ad("data/" + param + ".h5ad")
    adata = sort_cells(adata, adata_orig)
    # find cells that exist in orig dataset and this downsampled dataset

    if "pbmc4k" in param:
        t2g = pd.read_csv("data/homo_t2g.txt", delimiter="\t", index_col=1, header=None)
    else:
        t2g = pd.read_csv("data/mus_t2g.txt", delimiter="\t", index_col=1, header=None)

    temp = []
    for i in range(adata.var.index.shape[0]):
        if type(t2g.loc[adata.var.index[i], 2]) == pd.pandas.core.series.Series:
            temp.append(t2g.loc[adata.var.index[i], 2][0])
        else:
            temp.append(t2g.loc[adata.var.index[i], 2])
    adata.var = adata.var.set_index(np.array(temp))
    adata.obs_names_make_unique()
    # adata.write(adata_write_str)
    return adata
    # os.system("rm -d -r data/output_a")
    # os.system("rm -d -r data/output_b")


# %%


def main():
    # print("derp")
    # if(snakemake.params["input"] == "adata"):
    #     adata = methods.run(outputs = snakemake.output)
    # else:
    #     adata = methods.run(True,outputs = snakemake.output)
    adata = run_sampling(0, snakemake.wildcards[0])
    adata = methods.apply_method("downsample", adata, inputs=snakemake.input)
    adata.var["highly_variable"] = True
    if "diffexp" in snakemake.wildcards[0]:
        import scanpy as sc
        import numpy as np
        import pandas as pd

        np.savetxt("data/" + snakemake.wildcards[0] + "_leiden.txt", adata.obs["leiden"], fmt="%s")
        sc.tl.rank_genes_groups(adata, "leiden")
        pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(200).to_csv("data/" + snakemake.wildcards[0] + "_deg.csv", index=False)

    adata.write(snakemake.output[0])
    # adata_single.write(snakemake.output[1])


if __name__ == "__main__":
    main()
