#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 17:31:55 2023

@author: sinant
"""
import numpy as np
import pickle


# Read list to memory
def read_list(filename):
    # for reading also binary mode is important
    with open(filename, "rb") as fp:
        n_list = pickle.load(fp)
        return n_list


# %%
def parse_diffexp_full(file):
    return np.median(np.array(file, dtype=int), axis=0)
    # for i in range(len(file)):


def main():
    # %%
    # rr = read_list("data/diffexp_simul_pbmc_full.pickle")
    # arr_diffexp_full = np.array(diffexp_full)
    # diffexp_full = np.median(arr_diffexp_full, axis=0)

    # dd = read_list("data/diffexp_pbmc_full.pickle")
    diffexp_full = parse_diffexp_full(read_list("data/diffexp_pbmc_full.pickle"))
    diffexp_orig = parse_diffexp_full(read_list("data/diffexp_pbmc_orig_x.pickle"))
    diffexp_full_neuro = parse_diffexp_full(read_list("data/diffexp-neuro_full.pickle"))
    diffexp_orig_neuro = parse_diffexp_full(read_list("data/diffexp_neuro_orig_x.pickle"))

    diffexp_full_simul = parse_diffexp_full(read_list("data/diffexp-simul_pbmc_full.pickle"))
    diffexp_orig_simul = parse_diffexp_full(read_list("data/diffexp_simul_pbmc_orig_x.pickle"))

    diffexp_full_simul_neuro = parse_diffexp_full(read_list("data/diffexp_simul_neuro_full.pickle"))
    diffexp_orig_simul_neuro = parse_diffexp_full(read_list("data/diffexp_simul_neuro_orig_x.pickle"))

    # %%

    np.savetxt("data/diffexp_pbmc_full.csv", diffexp_full, delimiter=",", fmt="%s")
    np.savetxt("data/diffexp_pbmc_orig.csv", diffexp_orig, delimiter=",", fmt="%s")

    np.savetxt("data/diffexp_neuro_full.csv", diffexp_full_neuro, delimiter=",", fmt="%s")
    np.savetxt("data/diffexp_neuro_orig.csv", diffexp_orig_neuro, delimiter=",", fmt="%s")

    np.savetxt("data/diffexp_pbmc_full_simul.csv", diffexp_full_simul, delimiter=",", fmt="%s")
    np.savetxt("data/diffexp_pbmc_orig_simul.csv", diffexp_orig_simul, delimiter=",", fmt="%s")

    np.savetxt("data/diffexp_neuro_full_simul.csv", diffexp_full_simul_neuro, delimiter=",", fmt="%s")
    np.savetxt("data/diffexp_neuro_orig_simul.csv", diffexp_orig_simul_neuro, delimiter=",", fmt="%s")


if __name__ == "__main__":
    main()
