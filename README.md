# AM_2024
Analaysis for the batch correction paper

## Basic workflow

Using snakemake to create data-all.txt one can create all the data required for the analysis in the paper.

The basic workflow is as follows we call a single base module to create the base Anndata for each dataset. The run_pp in methods.py function is then called to preprocess each dataset. After this the Anndata for each particular method is created in R or python. The data is then corrected in an R script file or in apply_method() in methods.py. A batch corrected Anndata object is then created. Finally for each iteration we create the output files for the comparison metrics in create_output.py. We store them as .pickle files as the calculated arrays can differ in size. batch_report.py outputs those files to csv format to plot with the R file plot_results.R

For the differential expression results we run a similar workflow except the data for comparison is created using mod_mast.R from the Anndata files for each method. The pickle files are created with create_output_diffexp.py. Csv files used for plotting are then created with diffexp_report.py. The plots are then created with the plot_results.R file.

The following software/packages was used in the workflow. All packages from conda-forge, bioconda or default repositories.

* python 3.8.6
* scanpy 1.8.2
* leidenalg 0.8.9
* bbknn 1.4.1
* seurat 4.1.0
* harmony 0.1
* mast 1.20.0
* scvi-tools 0.15.0
* pyliger 0.2.0
* bustools 0.43.0
* kallisto 0.50.0

## Data

#### PBMC3K

3k PBMCs from a Healthy Donor - https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k

#### Mouse brain

1.3 Million Brain Cells from E18 Mice - https://www.10xgenomics.com/datasets/1-3-million-brain-cells-from-e-18-mice-2-standard-1-3-0

#### Mouse heart

1k Heart Cells from an E18 mouse (v3 chemistry) - https://www.10xgenomics.com/datasets/1-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0

#### PBMC4K

 4k PBMCs from a Healthy Donor - https://www.10xgenomics.com/datasets/4-k-pbm-cs-from-a-healthy-donor-2-standard-2-1-0

