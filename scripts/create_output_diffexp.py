#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import scanpy as sc
# import pandas as pd
import numpy as np

# from scipy.stats import mannwhitneyu
# from scipy.stats import ttest_ind

# from sklearn.metrics import pairwise_distances
# from scipy.spatial.distance import cdist
from os.path import exists
import pickle

# %%


def write_list(ls, filename):
    # store list in binary file so 'wb' mode
    file = []
    if exists(filename):
        # "with" statements are very handy for opening files.
        with open(filename, "rb") as fp:
            file = pickle.load(fp)
    file.append(ls)
    with open(filename, "wb") as fp:
        pickle.dump(file, fp)
        # print('Done writing list into a binary file')


# Read list to memory
def read_list(filename):
    # for reading also binary mode is important
    with open(filename, "rb") as fp:
        n_list = pickle.load(fp)
        return n_list


# %%


def main(inputs, wildcards, outputs):
    # compare total gene counts between B cells and cd8a t cells
    if "diffexp-pbmc" in wildcards[0]:
        p_val_string = "data/diffexp_pbmc_full.pickle"
        p_val_use_orig_x_string = "data/diffexp_pbmc_orig_x.pickle"
        # p_val_mast_string = "data/diffexp_pbmc_mast.pickle"
    if "diffexp-neuro" in wildcards[0]:
        p_val_string = "data/diffexp_neuro_full.pickle"
        p_val_use_orig_x_string = "data/diffexp_neuro_orig_x.pickle"
    if "simul-pbmc" in wildcards[0]:
        p_val_string = "data/diffexp_simul_pbmc_full.pickle"
        p_val_use_orig_x_string = "data/diffexp_simul_pbmc_orig_x.pickle"
        # p_val_mast_string = "data/diffexp_pbmc_mast.pickle"
    if "simul-neuro" in wildcards[0]:
        p_val_string = "data/diffexp_simul_neuro_full.pickle"
        p_val_use_orig_x_string = "data/diffexp_simul_neuro_orig_x.pickle"
        # p_val_mast_string = "data/diffexp_pbmc_mast.pickle"

    # read p files 1 and p files 2
    p_file_full = np.loadtxt(inputs[0], delimiter=" ")
    p_file_use_orig_clust = np.loadtxt(inputs[1], delimiter=" ")

    print("writing diffexp to pickle")
    write_list(p_file_full, p_val_string)
    write_list(p_file_use_orig_clust, p_val_use_orig_x_string)


if __name__ == "__main__":
    main(snakemake.input, snakemake.wildcards, snakemake.output)
    with open(snakemake.output[0], "w") as f:
        f.write("Create a new text file!")
