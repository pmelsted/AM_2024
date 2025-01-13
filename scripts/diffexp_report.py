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
    for data_str in snakemake.input:
        np.savetxt(f"{data_str[:-7]}.csv", parse_diffexp_full(read_list(data_str)), delimiter=",", fmt="%s")

    # %%



if __name__ == "__main__":
    main()
