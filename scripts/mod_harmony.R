library(harmony)
library(Seurat)

run <- function(){
  countsData <- read.table(file = snakemake@input[[1]], header = T, row.names=1, sep=",", as.is=T)
  pcaData <-  read.table(file = snakemake@input[[2]], header = T, row.names=1, sep=",", as.is=T)
  colnames(pcaData) <- paste0("PC_", 1:length(colnames(pcaData)))


  data <- CreateSeuratObject(counts = t(countsData), project = "data")
  data[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pcaData), key = "PC_", assay = DefaultAssay(data))

  
  batch <- read.csv(snakemake@input[[3]],header  = FALSE)
  Cells = data@meta.data
  Cells["batch"]= batch
  data <- AddMetaData(data,subset(Cells, select = c("batch")))
  data[["RNA"]]@scale.data = as.matrix(data[["RNA"]]@data)
  #variable_genes <- read.csv("data/variable_genes.csv", header = FALSE)
  #pbmc[["RNA"]]@var.features = variable_genes$V1
  #pbmc[["RNA"]]@scale.data = as.matrix(pbmc[["RNA"]]@data)
  
  data_res <- RunHarmony(data,"batch", plot_convergence = FALSE)
  res <- data_res[["harmony"]]@cell.embeddings
  write.table(t(as.matrix(res)),file=snakemake@output[[1]],row.names=FALSE,col.names=FALSE)
}

run()
