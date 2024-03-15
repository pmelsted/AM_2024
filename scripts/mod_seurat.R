library(Seurat)


run <- function(norm){
  if(norm==TRUE){
    countsData <- read.table(file = snakemake@input[[1]], header = T, row.names=1, sep=",", as.is=T)
    
    # ifnb.list <- lapply(X = pbmc.list, FUN = function(x) {
    #   x <- NormalizeData(x)
    # })
  }
  else{
    countsData <- read.table(file = snakemake@input[[1]], header = T, row.names=1, sep=",", as.is=T)

    
    #pbmc.list <- SplitObject(pbmc, split.by = "batch")
  }
  data <- CreateSeuratObject(counts = t(countsData), project = "data")
  batch <- read.csv(snakemake@input[[2]], header  = FALSE)
  Cells = data@meta.data
  Cells["batch"]= batch
  data <- AddMetaData(data,subset(Cells, select = c("batch")))
  
  #pbmc[["RNA"]]@var.features = variable_genes$V1
  data[["RNA"]]@scale.data = as.matrix(data[["RNA"]]@data)
  #pbmc[["RNA"]]@scale.data = as.matrix(pbmc[["RNA"]]@data)
  
  data.list <- SplitObject(data, split.by = "batch")
  if(norm==TRUE){
    variable_genes <- read.csv(snakemake@input[[3]], header = FALSE)$V1
    data.list <- lapply(X = data.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    k = gsub('-', '\\.', variable_genes)
    data.list$`0`@assays$RNA@var.features = k
    data.list$`1`@assays$RNA@var.features = k
    #features <- SelectIntegrationFeatures(object.list = pbmc.list)
    data.anchors <- FindIntegrationAnchors(data.list,scale=FALSE)
  }
  else{
    variable_genes <- read.csv(snakemake@input[[3]], header = FALSE)
    data[["RNA"]]@var.features = variable_genes$V1
    data.anchors <- FindIntegrationAnchors(data.list,scale=FALSE)#,anchor.features=variable_genes$V1)#,anchor.features = variable_genes)

  }
  
  
  
  # this command creates an 'integrated' data assay
  data.combined <- IntegrateData(anchorset = data.anchors)
  res <- data.combined[["integrated"]]@data
  #res = res[,order(colnames(res))]
  res = res[,match(colnames(data),colnames(res))]
  
  #tt =  t(as.matrix(res)) 
  #write.table(t(as.matrix(res)),file="data/r_to_adata_seurat.csv",row.names=FALSE,col.names=FALSE)
  if(norm==TRUE){
    write.table(t(as.matrix(res)),file=snakemake@output[[1]],row.names=TRUE,col.names=TRUE)
  } 
  else{
    write.table(t(as.matrix(res)),file=snakemake@output[[1]],row.names=TRUE,col.names=TRUE)
  }
    
  
}


if(as.logical(snakemake@params[["mode"]] == TRUE)){
  run(norm=TRUE)
}else{
  run(norm=FALSE)
}




