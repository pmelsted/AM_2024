#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
import pickle
import numpy as np
import scanpy as sc
import pandas as pd




# %matplotlib inline
# %%
methods_emb = ["combat","combatseq", "harmony", "liger","ligerv2" ,"mnn", "scvi", "seurat", "seuratv2"]
#methods_neuro_emb = ["combat","combatseq", "harmony", "liger", "mnn", "scvi", "seurat"]
methods_d_r = ["downsample"]#, "resample"]
methods = ["combat", "combatseq", "mnn", "scvi", "seurat", "seuratv2"]
methods_neuro = ["combat","combatseq" "mnn", "scvi", "seurat"]
#data_sets = ["pbmc3k","neuro","heart","pbmc4k","jejunum","simul_pbmc","simul_neuro"]
data_sets = ["pbmc3k","neuro","heart","pbmc4k","jejunum"]
# Read list to memory
def read_list(filename):
    """
    Read a list from a file.

    Args:
        filename (str): The path to the file.

    Returns:
        list: The list read from the file.
    """
    # for reading also binary mode is important
    with open(filename, "rb") as fp:
        n_list = pickle.load(fp)
        return n_list


def parse_and_plot_cc(data):
    """
    Parses the given data for the consesus clusters

    Parameters:
    - data: A list of lists representing the data.

    Returns:
    - The mean of the consensus clusters along the columns of the data.
    """
    arr = np.reshape(np.array(data), (len(data), len(data[0])))
    return np.mean(arr, axis=0)


def create_dict(data,emb=False, **kwargs):
    """
    Create a dictionary from the given data.
    Parameters:
    - data (list): A list of data.
    - emb (bool): Flag indicating whether the embedding is used or not. Default is False.
    - **kwargs: Additional keyword arguments.
    Returns:
    - data_melt (DataFrame): A melted DataFrame containing the transformed data.
    Raises:
    - None
    Example usage:
    data = [...] # list of data
    result = create_dict(data, emb=True, pbmc4k=True)
    """
    if not emb:

        if "pbmc4k" in kwargs or "heart" in kwargs:
            file_str = "pbmc4k" if "pbmc4k" in kwargs else "heart"
            adata = sc.read_h5ad(f"data/{file_str}-0_down.h5ad")
            in_all = adata.obs.index.str[:-2].values
            index = []
            for i in range(1, 25):
                adata = sc.read_h5ad(f"data/{file_str}-{i}_down.h5ad")
                in_all = np.intersect1d(in_all, adata.obs.index.str[:-2].values)

            for j in range(0, 25):
                adata = sc.read_h5ad(f"data/{file_str}-{j}_down.h5ad")
                index.append(np.in1d(adata.obs.index.str[:-2].values, in_all).nonzero()[0])
            index = np.array(index)

            down_rank = []
            for i in range(len(index)):
                down_rank.append(np.array(data[i][-2])[index[i].T])
            down_rank = np.median(np.array(down_rank), axis=0)

            # adata = sc.read_h5ad("data/{file_str}-0_resample.h5ad")
            # in_all = adata.obs.index.values
            # index = []
            # for i in range(1, 25):
            #     adata = sc.read_h5ad(f"data/{file_str}-{i}_resample.h5ad")
            #     in_all = np.intersect1d(in_all, adata.obs.index.values)
            # for j in range(0, 25):
            #     adata = sc.read_h5ad(f"data/{file_str}-{j}_resample.h5ad")
            #     index.append(np.in1d(adata.obs.index.values, in_all).nonzero()[0])
            # index = np.array(index)

            # resample_rank = []
            # for i in range(index.shape[0]):
            #     resample_rank.append(np.array(data[i][-1])[index[i].T])
            # resample_rank = np.median(np.array(resample_rank), axis=0)
            filt_data = [x[:len(methods)] for x in data]
            dic = {}
            for i, method in enumerate(methods):
                dic[method] = np.median(np.array([x[i] for x in filt_data]), axis=0)
            dic["downsample"] = down_rank
            #dic["resample"] = resample_rank

            #print(dic)
            df = pd.DataFrame({key: pd.Series(value) for key, value in dic.items()})
            #print(df)

            data_melt = pd.melt(
                df.reset_index(), id_vars=["index"], value_vars= methods + methods_d_r
            )
            #print(data_melt)

            data_melt = data_melt.dropna()
            return data_melt
    else:
        if "heart" in kwargs or "pbmc4k" in kwargs:
            file_str = "pbmc4k" if "pbmc4k" in kwargs else "heart"
            adata = sc.read_h5ad(f"data/{file_str}-0_down.h5ad")
            in_all = adata.obs.index.str[:-2].values
            index = []
            for i in range(1, 25):
                adata = sc.read_h5ad(f"data/{file_str}-{i}_down.h5ad")
                in_all = np.intersect1d(in_all, adata.obs.index.str[:-2].values)

            for j in range(0, 25):
                adata = sc.read_h5ad(f"data/{file_str}-{j}_down.h5ad")
                index.append(np.in1d(adata.obs.index.str[:-2].values, in_all).nonzero()[0])
            index = np.array(index)

            down_rank = []
            for i, val in enumerate(index):
                down_rank.append(np.array(data[i][-2])[val.T])
            down_rank = np.median(np.array(down_rank), axis=0)

            #adata = sc.read_h5ad(f"data/{file_str}-0_resample.h5ad")
            #in_all = adata.obs.index.values
            # index = []
            # for i in range(1, 25):
            #     adata = sc.read_h5ad(f"data/{file_str}-{i}_resample.h5ad")
            #     in_all = np.intersect1d(in_all, adata.obs.index.values)
            # for j in range(0, 25):
            #     adata = sc.read_h5ad(f"data/{file_str}-{j}_resample.h5ad")
            #     index.append(np.in1d(adata.obs.index.values, in_all).nonzero()[0])
            # index = np.array(index)

            # resample_rank = []
            # for i in range(index.shape[0]):
            #     resample_rank.append(np.array(data[i][-1])[index[i].T])
            # resample_rank = np.median(np.array(resample_rank), axis=0)
            dic = {}
            filt_data = [x[:len(methods_emb)] for x in data]
            for i, method in enumerate(methods_emb):
                dic[method] = np.median(np.array([x[i] for x in filt_data]), axis=0)
            dic["downsample"] = down_rank
            #dic["resample"] = resample_rank
            # print(dic)
            df = pd.DataFrame({key: pd.Series(value) for key, value in dic.items()})
            # print(df)

            data_melt = pd.melt(
                df.reset_index(),
                id_vars=["index"], 
                value_vars=methods_emb + methods_d_r,
            )
            # print(data_melt)
            data_melt = data_melt.dropna()
            return data_melt

def parse_rank(data, neuro=False, pbmc4k=False, heart=False):
    """
    Parses the given data and returns a melted DataFrame of nearest neighbor rank for x iterations.

    Parameters:
        data (list): The data to be parsed.
        neuro (bool, optional): If True, the data is parsed for neuro methods. Defaults to False.
        pbmc4k (bool, optional): If True, the data is parsed for pbmc4k methods. Defaults to False.
        heart (bool, optional): If True, the data is parsed for heart methods. Defaults to False.

    Returns:
        pandas.DataFrame: The melted DataFrame containing the parsed data of NN ranks, 
        median over x iterations.
    """
    if pbmc4k:
        data_melt = create_dict(data,emb=False,pbmc4k=True)
        return data_melt
    elif heart:
        data_melt = create_dict(data,emb=False, heart=True)
        return data_melt
    else:
        arr = np.array(data)
        arr = np.median(arr, axis=0)

        order = np.arange(0, np.array(data[0]).shape[1], 1)

        data_melt = {}
        if neuro:
            for i,method in enumerate(methods_neuro):
                data_melt[method] = arr[i]
            data_melt["order"] = order
            data_melt = pd.DataFrame                                                    (data_melt)
            data_melt = pd.melt(data_melt, id_vars=["order"], value_vars=methods_neuro)
        else:
            for i,method in enumerate(methods):
                data_melt[method] = arr[i]
            data_melt["order"] = order
            data_melt = pd.DataFrame(data_melt)
            data_melt = pd.melt(data_melt, id_vars=["order"], value_vars=methods)

        return data_melt


def parse_rank_emb(data,pbmc4k=False, heart=False,file_str=None):
    """
    Parses the given data and returns a melted DataFrame of nearest neighbor rank for the low dim embeddings.

    Parameters:
        data (list): The data to be parsed.
        neuro (bool, optional): If True, the data is parsed for neuro methods. Defaults to False.
        pbmc4k (bool, optional): If True, the data is parsed for pbmc4k methods. Defaults to False.
        heart (bool, optional): If True, the data is parsed for heart methods. Defaults to False.

    Returns:
        pandas.DataFrame: The melted DataFrame containing the parsed data of NN ranks for the low dim embeddings.
    """
    if pbmc4k:
        # print(data_melt)
        data_melt = create_dict(data,emb=True,pbmc4k=True)
        return data_melt
    if heart:
        data_melt = create_dict(data,emb=True,heart=True)
        # print(data_melt)
        return data_melt
    else:
        arr = np.array(data)
        arr = np.median(arr, axis=0)
        # clusters = []
        order = np.arange(0, np.array(data[0]).shape[1], 1)
        if False:
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
            data_melt = {}
            
            
            if "pbmc3k" not in file_str:
                if "ligerv2" in methods_emb:
                    methods_emb.remove("ligerv2")
            for i,method in enumerate(methods_emb):
                #print(i," ",method)
                data_melt[method] = arr[i]
            data_melt["order"] = order

            data_melt = pd.DataFrame(data_melt)
            data_melt = pd.melt(data_melt, id_vars=["order"], value_vars=methods_emb)

        return data_melt

    # %%


# %%


def main():
    cc_clust = []
    ari = []
    #for i in data_sets:
    for i in snakemake.input:
        pbmc4k = False
        heart = False
        if "emb" not in i:
            file_str = i[18:-7]
            #file_str = i
            if file_str == "pbmc4k":
                pbmc4k = True
            if file_str == "heart":
                heart = True
            print(file_str)
            parse_rank(read_list(
                f"data/nn_rank_file_{file_str}.pickle"),pbmc4k=pbmc4k,heart=heart).to_csv(f"data/{file_str}_plot.csv")
            # print(parse_and_plot_cc(read_list(f"data/cc_file_{file_str}.pickle")))
            cc_clust.append(parse_and_plot_cc(read_list(f"data/cc_file_{file_str}.pickle")))
            ari.append(parse_and_plot_cc(read_list(f"data/ari_{file_str}.pickle")))
            
        else:
            emb = True
            file_str = i[22:-7]
            if file_str == "pbmc4k":
                pbmc4k = True
            if file_str == "heart":
                heart = True
            print(file_str + " emb")
            parse_rank_emb(read_list(
                f"data/nn_rank_emb_file_{file_str}.pickle"),pbmc4k=pbmc4k,heart=heart,file_str = file_str).to_csv(f"data/{file_str}_emb_plot.csv")
    
        
    pd.DataFrame(cc_clust).to_csv("data/cc_all.csv")
    pd.DataFrame(ari).to_csv("data/ari_all.csv")
        


        
if __name__ == "__main__":
    main()
