#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import random

import methods

# import os
import string
import scanpy as sc
from copy import deepcopy

# %%


def id_generator(size=10, chars="ACGT"):
    return "".join(random.choice(chars) for _ in range(size))


# %%
def to_dic_np(ls):
    dic = {}
    dfa = ls.values
    # cells = np.array(dfa.index, type=object)
    for i in range(dfa.shape[0]):
        row = dfa[i, :]
        if row[0] in dic:
            if row[1] in dic[row[0]]:
                dic[row[0]][row[1]].append([row[2], row[3]])
            else:
                dic[row[0]][row[1]] = [[row[2], row[3]]]
        else:
            # dic[row["barcode"]] = [[row["umi"],row["ec"],row["count"]]]
            dic[row[0]] = {row[1]: [[row[2], row[3]]]}
    return dic


def to_dic(ls):
    dic = {}
    # dfa = ls.values
    for index, row in ls.iterrows():
        if row["barcode"] in dic:
            if row["umi"] in dic[row["barcode"]]:
                dic[row["barcode"]][row["umi"]].append([row["ec"], row["count"]])
            else:
                dic[row["barcode"]][row["umi"]] = [[row["ec"], row["count"]]]
        else:
            # dic[row["barcode"]] = [[row["umi"],row["ec"],row["count"]]]
            dic[row["barcode"]] = {row["umi"]: [[row["ec"], row["count"]]]}
    return dic


# %%


# %%


def resample_records(dic, cells):
    for barcode in dic.copy():
        if barcode in cells:
            for umi in dic[barcode].copy():
                for eq in dic[barcode][umi]:
                    c = np.random.default_rng().poisson(1, 1)[0]
                    if c == 0:
                        # do nothing
                        pass
                    else:
                        new_records = c * [eq]
                        for i in new_records:
                            new_umi = id_generator()
                            while new_umi in dic[barcode]:
                                new_umi = id_generator()
                            dic[barcode][new_umi] = [i]

    return dic


def to_bus_txt(dic):
    ls = []
    for barcode in dic:
        for umi in dic[barcode]:
            for eq in dic[barcode][umi]:
                ls.append([barcode, umi, eq[0], eq[1]])
    return pd.DataFrame(ls)


# %%


def run_sampling(param):
    ls = []
    res = []

    if "heart" in param:
        idx_file = "data/mus_index.idx"
        t2g_file = "data/mus_t2g.txt"
        prerun_cells = "data/heart_cells.txt"
        output_str = "data/resampled_heart/"
        adata_write_str = "data/adata_resampled_heart.h5ad"
        bus_out = "data/resampled_heart/output_sorted_corrected_heart.bus"
        txt_out = "data/resampled_heart/output_sorted_corrected_heart.txt"
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
            "cp -n mus_musculus/transcriptome.idx data/mus_index.idx",
            "cp -n mus_musculus/transcripts_to_genes.txt data/mus_t2g.txt"
            # "cd .. && rm -d mus_musculus.tar.gz"
        ]

        commands1 = [
            # kallisto_cmd,
            "bustools sort -o " + output_str + "output.s.bus " + output_str + "output.bus",
            "bustools correct -o " + output_str + "output.s.c.bus -w ./data/3M-february-2018.txt " + output_str + "output.s.bus",
            "bustools sort -o " + bus_out + " " + output_str + "output.s.c.bus"
            # bustools text -o ./data/output_sorted.txt " + bus_out
        ]
        if not os.path.isfile("data/3M-february-2018.txt"):
            os.system(
                "wget -O data/10x_version2_whitelist.txt https://github.com/BUStools/getting_started/releases/download/getting_started/10xv2_whitelist.txt"
            )

    elif "pbmc" in param:
        idx_file = "data/homo_index.idx"
        t2g_file = "data/homo_t2g.txt"
        output_str = "data/resampled_pbmc/"
        prerun_cells = "data/pbmc_cells.txt"
        adata_write_str = "data/adata_resampled_pbmc.h5ad"
        bus_out = "data/resampled_pbmc/output_sorted_corrected_pbmc.bus"
        txt_out = "data/resampled_pbmc/output_sorted_corrected_pbmc.txt"
        kallisto_cmd = "kallisto bus -o ./data/resampled_pbmc -i data/homo_index.idx -x 10xv2 -t 20 \
            ./data/fastqs/pbmc4k_S1_L001_R1_001.fastq.gz ./data/fastqs/pbmc4k_S1_L001_R2_001.fastq.gz  \
            ./data/fastqs/pbmc4k_S1_L002_R1_001.fastq.gz  ./data/fastqs/pbmc4k_S1_L002_R2_001.fastq.gz"
        commands0 = [
            "wget -nc https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz",
            "tar --skip-old-files -xf homo_sapiens.tar.gz",
            # "cd mus_musculus",
            "cp -n homo_sapiens/transcriptome.idx data/homo_index.idx",
            "cp -n homo_sapiens/transcripts_to_genes.txt data/homo_t2g.txt"
            # "cd .. && rm -d mus_musculus.tar.gz"
        ]

        commands1 = [
            # kallisto_cmd,
            "bustools sort -o data/resampled_pbmc/output.s.bus data/resampled_pbmc/output.bus",
            "bustools correct -o data/resampled_pbmc/output.s.c.bus -w ./data/10x_version2_whitelist.txt data/resampled_pbmc/output.s.bus",
            "bustools sort -o " + bus_out + " data/resampled_pbmc/output.s.c.bus"
            # bustools text -o ./data/output_sorted.txt " + bus_out
        ]
        if not os.path.isfile("data/10x_version2_whitelist.txt"):
            os.system(
                "wget -O data/10x_version2_whitelist.txt https://github.com/BUStools/getting_started/releases/download/getting_started/10xv2_whitelist.txt"
            )

    elif param == "mouse_neural":
        pass

    if not (os.path.isfile(idx_file) or os.path.isfile(t2g_file)):
        print("getting idx or t2g file")
        for c in commands0:
            os.system(c)

    if not os.path.isfile(output_str + "output.bus"):
        print("running kallisto")
        os.system(kallisto_cmd)

    if not os.path.isfile(bus_out):
        print("running sort correct sort")
        for c in commands1:
            os.system(c)
    if not os.path.isfile(txt_out):
        print("generating ")
        os.system("bustools text -o " + txt_out + " " + bus_out)

    with open(txt_out) as f:
        # gs = set()
        for line in f:
            # l = line.split()
            barcode, umi, ec, count = line.split()
            ec = int(ec)
            ls.append([barcode, umi, ec, count])

    # split data into 2 random equal size parts
    data = pd.read_csv(prerun_cells)
    data = data["barcode"]
    # batch_index = np.sort(np.random.choice(data.shape[0], size=int(0.5 * data.shape[0]), replace=False))
    batch = np.loadtxt("data/" + param + "_batch.csv").astype(int)

    batch_index = np.where(batch == 0)

    df_a = data.loc[batch_index]
    mask = np.ones(data.shape[0], dtype=bool)
    mask[batch_index] = False
    df_b = data[mask]

    reads = pd.DataFrame(ls, columns=["barcode", "umi", "ec", "count"])
    read_dic = to_dic_np(reads)
    print("starting resampling")
    read_dic_resample = resample_records(deepcopy(read_dic), np.array(df_a))

    to_txt = to_bus_txt(read_dic_resample)
    to_txt.columns = ["barcode", "umi", "ec", "count"]
    to_txt.to_csv(
        output_str + "output.resample.txt",
        sep="\t",
        columns=["barcode", "umi", "ec", "count"],
        header=False,
        index=False,
    )

    commands2 = [
        "bustools fromtext -o " + output_str + "resampled.bus " + output_str + "output.resample.txt",
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
        + "resampled.bus",
    ]

    # exec(open("scripts/sampling_bus.py").read())
    for c in commands2:
        print("reading in resampled and creating matrices")
        os.system(c)

    adata = sc.read(output_str + "output.mtx", cache=False)  # transpose the data*
    adata.var_names = pd.read_csv(
        output_str + "output.genes.txt",
        header=None,
    )[0]
    adata.obs_names = pd.read_csv(output_str + "output.barcodes.txt", header=None)[0]

    adata_orig = sc.read_h5ad("data/" + param + ".h5ad")
    same = np.intersect1d(adata_orig.obs.index.values, np.array(adata.obs.index.values))
    ind = np.in1d(adata_orig.obs.index.values, same).nonzero()[0]
    adata = adata[same, :]
    adata.obs["batch"] = batch[ind]
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

    return adata
    # os.system("rm -d -r data/output_a")
    # os.system("rm -d -r data/output_b")

    # pd.DataFrame(data=adata.layers["counts"].todense(), index=adata.obs_names, columns=adata.var_names).to_csv(outputs[3])
    # pd.DataFrame(adata.obs.index)
    # pd.DataFrame(adata.obs.index).to_csv("data/cells.txt",index=False)


def main():
    # print("derp")
    # if(snakemake.params["input"] == "adata"):
    #     adata = methods.run(outputs = snakemake.output)
    # else:
    #     adata = methods.run(True,outputs = snakemake.output)
    adata = run_sampling(snakemake.wildcards[0])
    # adata = run_sampling(0, snakemake.wildcards[0])
    adata = methods.apply_method("resample", adata, inputs=snakemake.input)
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
