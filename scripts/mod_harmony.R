library(harmony)
library(Seurat)

run <- function(){
  countsData <- read.table(file = snakemake@input[[1]], header = T, row.names=1, sep=",", as.is=T)
  pcaData <-  read.table(file = snakemake@input[[2]], header = T, row.names=1, sep=",", as.is=T)
  colnames(pcaData) <- paste0("PC_", 1:length(colnames(pcaData)))


  #emb = as.numeric(unlist(regmatches(snakemake@input[[1]], gregexpr("[[:digit:]]+", snakemake@input[[1]],))))[2]+1
  #nr = c(2,4,6,8,10,12,14,18,20,22,24,26,28,30,32,34,40,43,46,48,50,60,70,80,90,100)
  # print(emb)
  #print(nr[emb])
  #if(nr[emb] > 50){
  # pcaData <- cbind(pcaData,pcaData)
  #}
  
  data <- CreateSeuratObject(counts = t(countsData), project = "data")
  #data[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pcaData[,1:nr[emb]]), key = "PC_", assay = DefaultAssay(data))
  data[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pcaData), key = "PC_", assay = DefaultAssay(data))
  
  batch <- read.csv(snakemake@input[[3]],header  = FALSE)
  Cells = data@meta.data 
  Cells["batch"]= batch
  data <- AddMetaData(data,subset(Cells, select = c("batch")))
  data[["RNA"]]@scale.data = as.matrix(data[["RNA"]]@data)
  #variable_genes <- read.csv("data/variable_genes.csv", header = FALSE)
  #pbmc[["RNA"]]@var.features = variable_genes$V1
  #pbmc[["RNA"]]@scale.data = as.matrix(pbmc[["RNA"]]@data)
  
  data_res <- RunHarmony(data,"batch", dims.use=,plot_convergence = FALSE)
  res <- data_res[["harmony"]]@cell.embeddings
  write.table(t(as.matrix(res)),file=snakemake@output[[1]],row.names=FALSE,col.names=FALSE)
}

run()
