
if (!require("BiocManager", quietly = TRUE))
  #install.packages()
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install("sva")
library(sva)

run <- function(){
  countsData <- read.table(file = snakemake@input[[4]], header = T, row.names=1, sep=",", as.is=T)
  #pcaData <-  read.table(file = snakemake@input[[2]], header = T, row.names=1, sep=",", as.is=T)
  #colnames(pcaData) <- paste0("PC_", 1:length(colnames(pcaData)))

  countsData = t(countsData)

  
  #data <- CreateSeuratObject(counts = t(countsData), project = "data")
  #data[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pcaData[,1:nr[emb]]), key = "PC_", assay = DefaultAssay(data))
  #data[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pcaData), key = "PC_", assay = DefaultAssay(data))
  
  batch <- read.csv(snakemake@input[[2]],header  = FALSE)$V1
  #Cells = data@meta.data 
  #Cells["batch"]= batch
  #data <- AddMetaData(data,subset(Cells, select = c("batch")))
  #data[["RNA"]]@scale.data = as.matrix(data[["RNA"]]@data)
  #variable_genes <- read.csv("data/variable_genes.csv", header = FALSE)
  #pbmc[["RNA"]]@var.features = variable_genes$V1
  #pbmc[["RNA"]]@scale.data = as.matrix(pbmc[["RNA"]]@data)
  
  res <- ComBat_seq(countsData, batch=batch, group=NULL)
  
  write.table(t(as.matrix(res)),file=snakemake@output[[1]],row.names=FALSE,col.names=FALSE)
}

run()
