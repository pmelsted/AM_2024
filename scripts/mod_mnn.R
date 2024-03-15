

#library(tidyverse)
#library(remotes)
#library(readr)
#remotes::install_github("mojaveazure/seurat-disk")
#BiocManager::install("GenomeInfoDbData")
#BiocManager::install("GenomeInfoDb")

#library(batchelor)
library(Seurat)
#library(SeuratDisk)
library(batchelor)

run <- function(){
  countsData <- read.table(file = snakemake@input[[1]], header = T, row.names=1, sep=",", as.is=T)
  data <- CreateSeuratObject(counts = t(countsData), project = "data")
  batch <- read.csv(snakemake@input[[2]], 
                    header  = FALSE)
  Cells = data@meta.data
  Cells["batch"]= batch
  data <- AddMetaData(data,subset(Cells, select = c("batch")))
  out <- mnnCorrect(data[["RNA"]]@data,batch=data@meta.data["batch"]$batch)
  data[["RNA"]]@scale.data = out@assays@data@listData$corrected
  
  res <- data[["RNA"]]@scale.data
  write.table(t(res),file=snakemake@output[[1]],row.names=FALSE,col.names=FALSE)
}

run()
