#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
import numpy as np
import pickle
import scanpy as sc
import pandas as pd
from scipy.spatial import distance_matrix

# from data.methods import get_top_nn
# from data.methods import get_batch_loc_of_top_nn
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import *


# %matplotlib inline
# %%
# Read list to memory
def read_list(filename):
    # for reading also binary mode is important
    with open(filename, "rb") as fp:
        n_list = pickle.load(fp)
        return n_list


def parse_and_plot_cc(data):
    arr = np.reshape(np.array(data), (len(data), len(data[0])))
    return np.mean(arr, axis=0)


def parse_cc_new(cc_data):
    per_methods = []
    # if type(cc_data_temp != list):
    #     cc_data_temp = [cc_data_temp]
    cc_data_arr = np.array(cc_data_neuro).reshape((len(cc_data_neuro), len(cc_data_neuro[0])))
    return np.median(cc_data_arr, axis=0)


def parse_and_plot_rank(data, neuro=False, single=False, pbmc4k=False, heart=False):
    if pbmc4k:
        adata = sc.read_h5ad("data/pbmc4k-0_down.h5ad")
        in_all = adata.obs.index.str[:-2].values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{i}_down.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.str[:-2].values)

        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{j}_down.h5ad")
            index.append(np.in1d(adata.obs.index.str[:-2].values, in_all).nonzero()[0])
        index = np.array(index)

        down_rank = []
        for i in range(len(index)):
            down_rank.append(np.array(data[i][-2])[index[i].T])
        down_rank = np.median(np.array(down_rank), axis=0)

        adata = sc.read_h5ad("data/pbmc4k-0_resample.h5ad")
        in_all = adata.obs.index.values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{i}_resample.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.values)
        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{j}_resample.h5ad")
            index.append(np.in1d(adata.obs.index.values, in_all).nonzero()[0])
        index = np.array(index)

        resample_rank = []
        for i in range(index.shape[0]):
            resample_rank.append(np.array(data[i][-1])[index[i].T])
        resample_rank = np.median(np.array(resample_rank), axis=0)

        filt_data = [x[:5] for x in data]

        dic = {
            "combat": np.median(np.array([x[0] for x in filt_data]), axis=0),
            "mnn": np.median(np.array([x[1] for x in filt_data]), axis=0),
            "scvi": np.median(np.array([x[2] for x in filt_data]), axis=0),
            "seurat": np.median(np.array([x[3] for x in filt_data]), axis=0),
            "seuratv2": np.median(np.array([x[4] for x in filt_data]), axis=0),
            "downsample": down_rank,
            "resample": resample_rank,
            # "iter":np.median(np.array([x[0] for x in derp]),axis=0),
        }
        print(dic)
        df = pd.DataFrame({key: pd.Series(value) for key, value in dic.items()})
        print(df)

        data_melt = pd.melt(
            df.reset_index(), id_vars=["index"], value_vars=["combat", "mnn", "scvi", "seurat", "seuratv2", "downsample", "resample"]
        )
        # print(data_melt)

        data_melt = data_melt.dropna()
        print(data_melt)
        return data_melt
    elif heart:
        adata = sc.read_h5ad("data/heart-0_down.h5ad")
        in_all = adata.obs.index.str[:-2].values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/heart-{i}_down.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.str[:-2].values)

        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/heart-{j}_down.h5ad")
            index.append(np.in1d(adata.obs.index.str[:-2].values, in_all).nonzero()[0])
        index = np.array(index)

        down_rank = []
        for i in range(len(index)):
            down_rank.append(np.array(data[i][-2])[index[i].T])
        down_rank = np.median(np.array(down_rank), axis=0)

        adata = sc.read_h5ad("data/heart-0_resample.h5ad")
        in_all = adata.obs.index.values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/heart-{i}_resample.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.values)
        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/heart-{j}_resample.h5ad")
            index.append(np.in1d(adata.obs.index.values, in_all).nonzero()[0])
        index = np.array(index)

        resample_rank = []
        for i in range(index.shape[0]):
            resample_rank.append(np.array(data[i][-1])[index[i].T])
        resample_rank = np.median(np.array(resample_rank), axis=0)

        filt_data = [x[:5] for x in data]

        dic = {
            "combat": np.median(np.array([x[0] for x in filt_data]), axis=0),
            "mnn": np.median(np.array([x[1] for x in filt_data]), axis=0),
            "scvi": np.median(np.array([x[2] for x in filt_data]), axis=0),
            "seurat": np.median(np.array([x[3] for x in filt_data]), axis=0),
            "seuratv2": np.median(np.array([x[4] for x in filt_data]), axis=0),
            "downsample": down_rank,
            "resample": resample_rank,
            # "iter":np.median(np.array([x[0] for x in derp]),axis=0),
        }
        print(dic)
        df = pd.DataFrame({key: pd.Series(value) for key, value in dic.items()})
        print(df)

        data_melt = pd.melt(
            df.reset_index(), id_vars=["index"], value_vars=["combat", "mnn", "scvi", "seurat", "seuratv2", "downsample", "resample"]
        )
        # print(data_melt)

        data_melt = data_melt.dropna()
        print(data_melt)
        return data_melt
    else:
        arr = np.array(data)
        arr = np.median(arr, axis=0)

        # arr = np.zeros((data[:,:,0].shape[0],data[:,:,0].shape[1]))
        # for i in range(arr.shape[0]):
        #     for j in range(arr.shape[1]):
        #         arr[i,j] = np.median(data[i,j,:])

        clusters = []

        order = np.arange(0, np.array(data[0]).shape[1], 1)

        # if neuro:
        #     clusters = np.genfromtxt("data/leiden_neuro.csv", delimiter=",", skip_header=1)[:, 1]

        # else:
        #     clusters = np.genfromtxt("data/leiden.csv", delimiter=",", skip_header=1)[:, 1]
        if neuro:
            data = {"combat": arr[0], "mnn": arr[1], "scvi": arr[2], "seurat": arr[3], "order": order}
            data = pd.DataFrame(data)
            data_melt = pd.melt(data, id_vars=["order"], value_vars=["combat", "mnn", "scvi", "seurat"])
        else:
            data = {"combat": arr[0], "mnn": arr[1], "scvi": arr[2], "seurat": arr[3], "seuratv2": arr[4], "order": order}
            data = pd.DataFrame(data)
            data_melt = pd.melt(data, id_vars=["order"], value_vars=["combat", "mnn", "scvi", "seurat", "seuratv2"])

        return data_melt


def parse_and_plot_rank_emb(data, neuro=False, pbmc4k=False, heart=False):
    if pbmc4k:
        adata = sc.read_h5ad("data/pbmc4k-0_down.h5ad")
        in_all = adata.obs.index.str[:-2].values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{i}_down.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.str[:-2].values)

        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{j}_down.h5ad")
            index.append(np.in1d(adata.obs.index.str[:-2].values, in_all).nonzero()[0])
        index = np.array(index)

        down_rank = []
        for i in range(len(index)):
            down_rank.append(np.array(data[i][-2])[index[i].T])
        down_rank = np.median(np.array(down_rank), axis=0)

        adata = sc.read_h5ad("data/pbmc4k-0_resample.h5ad")
        in_all = adata.obs.index.values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{i}_resample.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.values)
        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/pbmc4k-{j}_resample.h5ad")
            index.append(np.in1d(adata.obs.index.values, in_all).nonzero()[0])
        index = np.array(index)

        resample_rank = []
        for i in range(index.shape[0]):
            resample_rank.append(np.array(data[i][-1])[index[i].T])
        resample_rank = np.median(np.array(resample_rank), axis=0)

        filt_data = [x[:7] for x in data]

        dic = {
            "combat": np.median(np.array([x[0] for x in filt_data]), axis=0),
            "harmony": np.median(np.array([x[1] for x in filt_data]), axis=0),
            "liger": np.median(np.array([x[2] for x in filt_data]), axis=0),
            "mnn": np.median(np.array([x[3] for x in filt_data]), axis=0),
            "scvi": np.median(np.array([x[4] for x in filt_data]), axis=0),
            "seurat": np.median(np.array([x[5] for x in filt_data]), axis=0),
            "seuratv2": np.median(np.array([x[6] for x in filt_data]), axis=0),
            "downsample": down_rank,
            "resample": resample_rank,
            # "iter":np.median(np.array([x[0] for x in derp]),axis=0),
        }
        # print(dic)
        df = pd.DataFrame({key: pd.Series(value) for key, value in dic.items()})
        # print(df)

        data_melt = pd.melt(
            df.reset_index(),
            id_vars=["index"],
            value_vars=["combat", "harmony", "liger", "mnn", "scvi", "seurat", "seuratv2", "downsample", "resample"],
        )
        # print(data_melt)

        data_melt = data_melt.dropna()
        # print(data_melt)
        return data_melt
    elif heart:
        adata = sc.read_h5ad("data/heart-0_down.h5ad")
        in_all = adata.obs.index.str[:-2].values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/heart-{i}_down.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.str[:-2].values)

        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/heart-{j}_down.h5ad")
            index.append(np.in1d(adata.obs.index.str[:-2].values, in_all).nonzero()[0])
        index = np.array(index)

        down_rank = []
        for i in range(len(index)):
            down_rank.append(np.array(data[i][-2])[index[i].T])
        down_rank = np.median(np.array(down_rank), axis=0)

        adata = sc.read_h5ad("data/heart-0_resample.h5ad")
        in_all = adata.obs.index.values
        index = []
        for i in range(1, 25):
            adata = sc.read_h5ad(f"data/heart-{i}_resample.h5ad")
            in_all = np.intersect1d(in_all, adata.obs.index.values)
        for j in range(0, 25):
            adata = sc.read_h5ad(f"data/heart-{j}_resample.h5ad")
            index.append(np.in1d(adata.obs.index.values, in_all).nonzero()[0])
        index = np.array(index)

        resample_rank = []
        for i in range(index.shape[0]):
            resample_rank.append(np.array(data[i][-1])[index[i].T])
        resample_rank = np.median(np.array(resample_rank), axis=0)

        filt_data = [x[:7] for x in data]
        print(len(filt_data))
        dic = {
            "combat": np.median(np.array([x[0] for x in filt_data]), axis=0),
            "harmony": np.median(np.array([x[1] for x in filt_data]), axis=0),
            "liger": np.median(np.array([x[2] for x in filt_data]), axis=0),
            "mnn": np.median(np.array([x[3] for x in filt_data]), axis=0),
            "scvi": np.median(np.array([x[4] for x in filt_data]), axis=0),
            "seurat": np.median(np.array([x[5] for x in filt_data]), axis=0),
            "seuratv2": np.median(np.array([x[6] for x in filt_data]), axis=0),
            "downsample": down_rank,
            "resample": resample_rank,
            # "iter":np.median(np.array([x[0] for x in derp]),axis=0),
        }
        # print(dic)
        df = pd.DataFrame({key: pd.Series(value) for key, value in dic.items()})
        # print(df)

        data_melt = pd.melt(
            df.reset_index(),
            id_vars=["index"],
            value_vars=["combat", "harmony", "liger", "mnn", "scvi", "seurat", "seuratv2", "downsample", "resample"],
        )
        # print(data_melt)

        data_melt = data_melt.dropna()
        # print(data_melt)
        return data_melt
    else:
        arr = np.array(data)
        arr = np.median(arr, axis=0)
        # clusters = []
        order = np.arange(0, np.array(data[0]).shape[1], 1)
        if arr.shape[0] == 8:
            data = {
                "combat": arr[0],
                "harmony": arr[1],
                "liger": arr[2],
                "ligerv2": arr[3],
                "mnn": arr[4],
                "scvi": arr[5],
                "seurat": arr[6],
                "seuratv2": arr[7],
                "order": order,
            }
            data = pd.DataFrame(data)
            data_melt = pd.melt(
                data, id_vars=["order"], value_vars=["combat", "harmony", "liger", "ligerv2", "mnn", "scvi", "seurat", "seuratv2"]
            )
        else:
            data = {
                "combat": arr[0],
                "harmony": arr[1],
                "liger": arr[2],
                "mnn": arr[3],
                "scvi": arr[4],
                "seurat": arr[5],
                "seuratv2": arr[6],
                "order": order,
            }
            data = pd.DataFrame(data)
            data_melt = pd.melt(data, id_vars=["order"], value_vars=["combat", "harmony", "liger", "mnn", "scvi", "seurat", "seuratv2"])

        return data_melt

    # %%


# %%


def main():
    for i in [
        read_list("data/cc_file.pickle"),
        read_list("data/cc_file_neuro.pickle"),
        read_list("data/cc_file_heart.pickle"),
        read_list("data/cc_file_pbmc4k.pickle"),
        read_list("data/cc_file_simul_pbmc.pickle"),
        read_list("data/cc_file_simul_neuro.pickle"),
    ]:
        print(parse_and_plot_cc(i))

    # pbmc3k TO csv
    parse_and_plot_rank(read_list("data/nn_rank_file.pickle")).to_csv("data/pbmc_plot.csv")
    parse_and_plot_rank_emb(read_list("data/nn_rank_emb_file.pickle")).to_csv("data/pbmc_emb_plot.csv")

    # neuro save to csv
    parse_and_plot_rank(read_list("data/nn_rank_file_neuro.pickle")).to_csv("data/neuro_plot.csv")
    parse_and_plot_rank_emb(read_list("data/nn_rank_emb_file_neuro.pickle")).to_csv("data/neuro_emb_plot.csv")

    # pbmc4k save to csv
    parse_and_plot_rank(read_list("data/nn_rank_file_pbmc4k.pickle"), pbmc4k=True).to_csv("data/pbmc4k_plot.csv")
    parse_and_plot_rank_emb(read_list("data/nn_rank_emb_file_pbmc4k.pickle"), pbmc4k=True).to_csv("data/pbmc4k_emb_plot.csv")

    # heart save to csv

    parse_and_plot_rank(read_list("data/nn_rank_file_heart.pickle"), heart=True).to_csv("data/heart_plot.csv")
    parse_and_plot_rank_emb(read_list("data/nn_rank_emb_file_heart.pickle"), heart=True).to_csv("data/heart_emb_plot.csv")

    # pbmc simulto csv
    parse_and_plot_rank(read_list("data/nn_rank_file_simul_pbmc.pickle"), resamp=True).to_csv("data/simul_pbmc_plot.csv")
    parse_and_plot_rank_emb(read_list("data/nn_rank_emb_file_simul_pbmc.pickle"), resamp=True).to_csv("data/simul_pbmc_emb_plot.csv")

    # pbmc simulto csv
    parse_and_plot_rank(read_list("data/nn_rank_file_simul_neuro.pickle"), resamp=True).to_csv("data/simul_neuro_plot.csv")
    parse_and_plot_rank_emb(read_list("data/nn_rank_emb_file_simul_neuro.pickle"), resamp=True).to_csv("data/simul_neuro_emb_plot.csv")


if __name__ == "__main__":
    main()
