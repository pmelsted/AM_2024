library(Seurat)
library(MAST)

run <- function(){
  adata <- create_data()
  adata_bbknn <- create_data("bbknn")
  adata_combat <- create_data("combat")
  adata_harmony <- create_data("harmony")
  adata_liger <- create_data("liger")
  adata_mnn <- create_data("mnn")
  adata_scvi <- create_data("scvi")
  adata_seurat <- create_data("seurat")
  
  
  methods <- list(adata,adata_bbknn,adata_combat,adata_harmony,adata_liger,adata_mnn,adata_scvi,adata_seurat)
  
  batch_clust_data <- test_batch_cl_data(methods)
  
  orig_clust_data <- test_orig_cl_data(list(adata,adata_bbknn,adata_combat,adata_harmony,adata_liger,adata_mnn,adata_scvi,adata_seurat),adata)
  
  write.table(batch_clust_data,file=snakemake@output[[1]],row.names=FALSE,col.names=FALSE)
  write.table(orig_clust_data,file=snakemake@output[[2]],row.names=FALSE,col.names=FALSE)
}
#workflow stages

#read in orig csv as in function above

#read in leiden clusters 

# read in deg files



create_data <- function(name=""){
  
  # get right X array, leiden clustering, deg for each cluster
  # for each methods
  batch_str <- paste("data/",snakemake@wildcards[[1]],"_batch.csv",sep = "")
  if(name == ""){
    file_str <- paste("data/",snakemake@wildcards[[1]],"_to_csv.csv",sep = "")
    leiden_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_leiden.txt",sep = "")
    diff_exp_gene_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_deg.csv",sep = "")
    
  }
  else if(name == "liger"){
    file_str <- paste("data/",snakemake@wildcards[[1]],"_to_csv.csv",sep = "")
    leiden_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_leiden.txt",sep = "")
    #diff_exp_gene_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_deg.csv",sep = "")
  }
  else if(name == "bbknn" | name == "harmony"){
    file_str <- paste("data/",snakemake@wildcards[[1]],"_to_csv.csv",sep = "")
    leiden_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_leiden.txt",sep = "")
    diff_exp_gene_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_deg.csv",sep = "")
  }
  else{
    file_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_X.csv",sep = "")
    leiden_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_leiden.txt",sep = "")
    diff_exp_gene_clusters_str <- paste("data/",snakemake@wildcards[[1]],"_",name,"_deg.csv",sep = "")
  }

  countsData <- read.table(file = file_str, header = T, row.names=1, sep=",", as.is=T)
  if(grepl("pbmc",snakemake@wildcards[[1]])){
    proj = "pbmc"
  }
  if(grepl("neuro",snakemake@wildcards[[1]])){
    proj = "neuro"
  }
  data <- CreateSeuratObject(counts = t(countsData), project = proj)
  
  batch <- read.csv(batch_str,header  = FALSE)
  Cells = data@meta.data
  Cells["batch"]= batch

  data <- AddMetaData(data,subset(Cells, select = c("batch")))
  data[["RNA"]]@scale.data = as.matrix(data[["RNA"]]@data)
  
  
  data.sce <- as.SingleCellExperiment(data)
  
  
  cl = read.csv(leiden_clusters_str,header=FALSE)$V1
  colData(data.sce)$ident <- as.factor(cl)




  # genes <- apply(X=genes,MARGIN=c(1,2),FUN=tolower)
  # clust_b <- which(genes == "ms4a1",arr.ind=T)[2]
  # clust_t <- which(genes == "cd8a",arr.ind=T)[2]
  #get cell in that clust
  
  data.sca <- SceToSingleCellAssay(data.sce, class = "SingleCellAssay", check_sanity = TRUE)
  #data.sca@int_metadata$clusters <- c(clust_b,clust_t)
  if(name == "liger"){
    d = read.csv(paste("data/",snakemake@wildcards[[1]],"_",name,"_","clust_pair.txt",sep = ""),header=FALSE)
    
    
    #d <- read.csv(paste("data/","liger_clust_pair.txt",sep = ""),header=FALSE)
    data.sca@int_metadata$clusters <- c(d$V1[1],d$V1[2])
  }
  else{
    genes <- read.csv(diff_exp_gene_clusters_str)
    data.sca@int_metadata$clusters <- get_genes(genes)
  }

  return(data.sca)
}

get_genes <- function(diff_exp_genes){
  genes <- apply(X=diff_exp_genes,MARGIN=c(1,2),FUN=tolower)
  if(grepl("pbmc",snakemake@wildcards[[1]])){
    gene1 = "ms4a1"
    gene2 = "cd8a"
  }
  if(grepl("neuro",snakemake@wildcards[[1]])){
    gene1 = "myl9"
    gene2 = "ctss"
  }
  gene_search1 <- which(genes == gene1,arr.ind=T)
  gene_search2 <- which(genes == gene2,arr.ind=T)
  
  if(length(gene_search1) == 2){
    clust_1 <- gene_search1[2]
  }
  else{
    clust_1 <- as.numeric(gene_search1[1,2])
  }
  if(length(gene_search2) == 2){
    clust_2 <- gene_search2[2]
  }
  else{
    clust_2 <- as.numeric(gene_search2[1,2])
  }
  return(c(clust_1-1,clust_2-1))
}
add_id_cl_col <- function(data,cl){
  colData(data)["ident_select"] <- with(colData(data), ifelse(ident==toString(cl),1,0))
  # colData(data)["ident_select"] <- data.frame(colData(data)) %>% 
  #   mutate("is_ident" = as.factor(if_else(ident==toString(cl),1,0))) %>% 
  #   select("is_ident")
  return(data)
}

get_p_vals <- function(data,batch=FALSE,cluster=FALSE){
  if(batch == TRUE & cluster == FALSE){
    #print("this is run")
    zlm.output <- zlm(~batch, data,)
    zlm.lr_clust <- lrTest(zlm.output, "batch")
  }
  if(batch == FALSE & cluster == TRUE){
    zlm.output <- zlm(~ident_select, data,)
    zlm.lr_clust <- lrTest(zlm.output, "ident_select")
  }
  if(batch & cluster){
    zlm.output <- zlm(~ident_select + batch, data,)
    zlm.lr_clust <- lrTest(zlm.output, CoefficientHypothesis(c('batch','ident_select')))
  }
  
  
  names <- rownames(zlm.lr_clust)
  values <- zlm.lr_clust[,3,3]
  res <- data.frame(values,row.names = names)
  return(res)
  
}

#filter cells and genes
#get p vals for LRT of batch coef model vs intercept only model
#get p vals for LRT of cluster coef model vs intercept only model
#get p vals for LRT of cluster coef in cluster + batch model vs intercept model.
# return how many genes were statistically significant after multiple test correction
# we use bonferroni because we are interested in cut off
# of family wide error rate not trying to find best

test_batch_cl_data <- function(methods){
  ind = 0
  gene_p_vals = NULL
  orig_genes = NULL
  for(data.sca in methods){
    
    
    clust_1 <- data.sca@int_metadata$clusters[1]
    clust_2 <- data.sca@int_metadata$clusters[2]
    #filter data to only use cells types of interest
    data.sca_filt <- data.sca[,colData(data.sca)$ident == toString(clust_1) | colData(data.sca)$ident == toString(clust_2) ]
    #filter genes used to identify cluster of b cells and t cells
    # so as to not bias model
    if(grepl("pbmc",snakemake@wildcards[[1]])){
      gene1 = "ms4a1"
      gene2 = "cd8a"
    }
    if(grepl("neuro",snakemake@wildcards[[1]])){
      gene1 = "myl9"
      gene2 = "ctss"
    }
    index1 <- which(tolower(rowData(data.sca)$primerid)== gene1)
    index2 <- which(tolower(rowData(data.sca)$primerid)== gene2)
    if(length(index1) == 1){
      data.sca_filt[-index1,]
    }
    if(length(index2) == 1){
      data.sca_filt[-index2,]
    }
    
    data.sca_filt <- add_id_cl_col(data.sca_filt,clust_1)
    
    p_vals1 <- get_p_vals(data.sca_filt,clust=T)
    orig_nr <- length(p_vals1[p_vals1 < 0.05/length(p_vals1$values)]) #179 mode <l with cluster coefficient, orig
  
    
    p_vals2 <- get_p_vals(data.sca_filt,batch=T,clust=F)
    batch_nr <- length(p_vals2[p_vals2 < 0.05/length(p_vals2$values)]) #0 model with batch coefficient, orig
    #get same info but corrected for batch only
    
    p_vals3 <- get_p_vals(data.sca_filt,clust=T,batch=T) # 176 model with batch + cluster
    batch_clust_nr <- length(p_vals3[p_vals3 < 0.05/length(p_vals3$values)]) 
    
    
    #ratio of same genes
    
    
    #match is done on model with only cluster as covariate
    
    if(ind == 0){
      gene_p_vals = data.frame(mth_0 = c(orig_nr,batch_nr,batch_clust_nr,0))#,row.names = rownames(p_vals))
      p_vals_df = data.frame(get_p_vals(data.sca_filt,clust=T))
      orig_genes = rownames(p_vals_df[p_vals_df$values < 0.05/length(p_vals_df$values),,drop=FALSE])
      #orig_genes = rownames(p_vals1[p_vals1 < 0.05/length(p_vals1$values)])
      #orig_genes = rownames(get_p_vals(data.sca_filt,clust=T))
    }
    else{
      p_vals_df = data.frame(p_vals1)
      other_method_genes = rownames(p_vals_df[p_vals_df$values < 0.05/length(p_vals_df$values),,drop=FALSE])
      #other_method_genes = rownames(p_vals1[p_vals1 < 0.05/length(p_vals1$values)])
      nr_match_genes <- length(intersect(orig_genes, other_method_genes))
      gene_p_vals[paste("mth_",toString(ind),sep = "")] <- c(orig_nr,batch_nr,batch_clust_nr,nr_match_genes)
      
      
      
    }
    
    ind = ind + 1
  }
  
  return(gene_p_vals)
  
}

test_orig_cl_data <- function(methods,orig){
  ind = 0
  gene_p_vals = NULL
  for(data.sca in methods){
    
    #the 2 cell clusters to compare
    clust_1 <- data.sca@int_metadata$clusters[1]
    clust_2 <- data.sca@int_metadata$clusters[2]
    #filter data to only use cells types of interest
    data.sca_filt <- orig
    colData(data.sca_filt)$ident <- colData(data.sca)$ident
    #colData(data.sca_filt)$ident <- colData(orig)$ident
    data.sca_filt <- data.sca_filt[,colData(data.sca_filt)$ident == toString(clust_1) | colData(data.sca_filt)$ident == toString(clust_2) ]
    #filter genes used to identify cluster of b cells and t cells
    # so as to not bias model
    if(grepl("pbmc",snakemake@wildcards[[1]])){
      gene1 = "ms4a1"
      gene2 = "cd8a"
    }
    if(grepl("neuro",snakemake@wildcards[[1]])){
      gene1 = "myl9"
      gene2 = "ctss"
    }
    index1 <- which(tolower(rowData(data.sca_filt)$primerid)== gene1)
    index2 <- which(tolower(rowData(data.sca_filt)$primerid)== gene2)
    if(length(index1) == 1){
      data.sca_filt[-index1,]
    }
    if(length(index2) == 1){
      data.sca_filt[-index2,]
    }
    
    data.sca_filt <- add_id_cl_col(data.sca_filt,clust_1)
    
    p_vals1 <- get_p_vals(data.sca_filt,clust=T)
    orig_nr <- length(p_vals1[p_vals1 < 0.05/length(p_vals1$values)]) #179 mode <l with cluster coefficient, orig
    
    p_vals2 <- get_p_vals(data.sca_filt,batch=T,clust=F)
    batch_nr <- length(p_vals2[p_vals2 < 0.05/length(p_vals2$values)]) #0 model with batch coefficient, orig
    #get same info but corrected for batch only
    
    p_vals3 <- get_p_vals(data.sca_filt,clust=T,batch=T) # 176 model with batch + cluster
    batch_clust_nr <- length(p_vals3[p_vals3 < 0.05/length(p_vals3$values)]) 
    
    if(ind == 0){
      gene_p_vals = data.frame(mth_0 = c(orig_nr,batch_nr,batch_clust_nr,0))#,row.names = rownames(p_vals))
      p_vals_df = data.frame(get_p_vals(data.sca_filt,clust=T))
      orig_genes = rownames(p_vals_df[p_vals_df$values < 0.05/length(p_vals_df$values),,drop=FALSE])
    }
    else{
      p_vals_df = data.frame(p_vals1)
      other_method_genes = rownames(p_vals_df[p_vals_df$values < 0.05/length(p_vals_df$values),,drop=FALSE])
      nr_match_genes <- length(intersect(orig_genes,other_method_genes))
      gene_p_vals[paste("mth_",toString(ind),sep = "")] <- c(orig_nr,batch_nr,batch_clust_nr,nr_match_genes)
    }
    
    ind = ind + 1
  }
  
  return(gene_p_vals)
  
}





run()
