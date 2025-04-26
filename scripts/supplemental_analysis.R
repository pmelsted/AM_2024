
#' ---
#'  
#'  
#'  title: "Supplemental Analysis"
#'  Run the Seurat preprocessing for the PBMC3K dataset
#'  Takes the raw values and performs the preprocessing steps
#'  Finding variable features and normalizing the data 
run_seurat_preprocessing <- function(raw_str,batch_str,count_str){
    data <- read.table(file = raw_str, header = T, row.names=1, sep=",", as.is=T)
    data <- CreateSeuratObject(counts = t(data), project = "data")
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    
    data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    
    batch <- read.csv(batch_str,header  = FALSE)
    Cells = data@meta.data
    Cells["batch"]= batch
    data <- AddMetaData(data,subset(Cells, select = c("batch")))
    
    data.list <- SplitObject(data, split.by = "batch")
    
    # normalize and identify variable features for each dataset independently
    data.list <- lapply(X = data.list, FUN = function(x) {
        x <- NormalizeData(x,verbose=FALSE)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = data.list,verbose=FALSE)
    
    data.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features,verbose=FALSE)
    data.combined <- IntegrateData(anchorset = data.anchors,verbose=FALSE)
    
    
    X = t(as.matrix(data.combined@assays$integrated@data))
    
    write.table(t(as.matrix(data.combined@assays$integrated@data)),file=paste(count_str,sep=""),row.names=TRUE,col.names=TRUE)
    

}

#theta change
#'  
#'  Run the Harmony batch correction and vary the theta parameter. 
#'  Get NN rank change and return the result.
harmony_theta <- function(pca_str){
    
    pca_nn  <-  read.table(file = pca_str, header = T, row.names=1, sep=",", as.is=T)
    pca_nn <- dist(pca_nn, method="euclidean")
    pca_nn <- as.matrix(pca_nn)
    pca_nn = get_top_k_neighbors(pca_nn, 2637)
    ans <- numeric(100)
    for (i in 1:100){
        tryCatch(
            {
                theta= 2*i
                sigma = 0.1
                lambda = 1
                X_harm <- run_harm(theta,sigma,lambda)
                dist_harm <- dist(X_harm, method="euclidean")
                dist_harm <- as.matrix(dist_harm)
                top_10_nn = get_top_k_neighbors(dist_harm, 30)
                res <- nn_rank_change(top_10_nn, pca_nn)
                rank_change_theta1 <-   pmap(as.data.frame(res),function(...) median(c(...)))
                val <- median(map_dbl(rank_change_theta1,pluck))
                ans[i] = val 
            },
            error = function(e) {
                print(e)
                return(ans)
            },
            finally={
                #return(val)
            }
        )

    }
    return(ans)
    
}


#theta change
#'  
#'  Run the Harmony batch correction and vary the sigma parameter. 
#'  Get NN rank change and return the result.
harmony_sigma <- function(){
    
    
    pca_nn  <-  read.table(file = pca_str, header = T, row.names=1, sep=",", as.is=T)
    pca_nn <- dist(pca_nn, method="euclidean")
    pca_nn <- as.matrix(pca_nn)
    pca_nn = get_top_k_neighbors(pca_nn, 2637)
    ans <- numeric(50)
    j = seq(from = 0.0001, to = 0.5, length.out = 50)
    for (i in 1:50){
        tryCatch(
            {
                theta= 0
                sigma = j[i]
                lambda = 0.0001
                X_harm <- run_harm(theta,sigma,lambda)
                dist_harm <- dist(X_harm, method="euclidean")
                dist_harm <- as.matrix(dist_harm)
                top_10_nn = get_top_k_neighbors(dist_harm, 30)
                res <- nn_rank_change(top_10_nn, pca_nn)
                rank_change_theta1 <-   pmap(as.data.frame(res),function(...) median(c(...)))
                val <- median(map_dbl(rank_change_theta1,pluck))
                ans[i] = val 
                print(val)
            },
            error = function(e) {
                print(e)
                return(ans)
            },
            finally={
                #return(val)
            }
        )
    }
    return(ans)
    
}

