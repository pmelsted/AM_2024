library(MAST)
library(sva)
library(Seurat)
library(tidyverse)
#VOLCANO
create_sca <- function(d){
  dat = read.table(file = d, header = T, row.names=1, sep=",", as.is=T)
  #print(substring(d,15,nchar(d)-4))
  cl = read.csv(paste("data/pbmc3k-0_",substring(d,15,nchar(d)-4),"_leiden.txt",sep=""),header=FALSE)$V1
  data <- CreateSeuratObject(counts = t(dat), project = "data")
  batch <- read.csv("data/pbmc3k-0_batch.csv", header  = FALSE)$V1
  clust = read.csv("data/leiden_pbmc3k.csv",header=FALSE,skip=1)$V2
  #data <- CreateSeuratObject(counts = t(countsData), project = proj)
  Cells = data@meta.data
  Cells["batch"]= batch
  Cells["leiden"] = as.factor(cl)
  #Cells["ident"] = as.factor(clust)
  data <- AddMetaData(data,subset(Cells, select = c("batch","leiden")))
  Idents(data) <- "leiden"
  #data <- AddMetaData(data,subset(Cells, select = c("batch")))
  #colData(data.sce)$ident <- as.factor(clust)
  #pbmc[["RNA"]]@var.features = variable_genes$V1
  data@misc <- c(data@misc,get_clusts(substring(d,15,nchar(d)-4)))
  data[["RNA"]]$data = data[["RNA"]]$counts
  #batch <- read.csv(batch_str,header  = FALSE)
  
  data.sce <- as.SingleCellExperiment(data)
  assay(data.sce,"data") <- counts(data.sce)
  assay(data.sce,"scale.data") <- counts(data.sce)
  
  #cl = read.csv(leiden_clusters_str,header=FALSE)$V1
  colData(data.sce)$orig.ident <- as.factor(clust)
  
  if(substring(d,15,nchar(d)-4) == "scvi"){
    logcounts(data.sce) <- log2(counts(data.sce) + 1)
  }
  
  data.sca <- SceToSingleCellAssay(data.sce, class = "SingleCellAssay", check_sanity = TRUE)
  data.sca@int_metadata$c = get_clusts(substring(d,15,nchar(d)-4))
  data.sca@int_metadata$name <- substring(d,15,nchar(d)-4)
  return(data.sca)
  
}
read_data_sca <- function(files,file_str) {
  # Create the full file paths based on the directory and specified file names
  file_names = lapply(files,function(x) paste(file_str,x,".csv",sep=""))
  file_paths <- file.path("data", file_names)
  
  # Read each specified file and store it in a list
  #data_list <- lapply(file_paths, function(x) )
  data_list <- lapply(file_paths,create_sca)
  #data_list <- lapply(data_list,function(x) x[["RNA"]]$data <- x[["RNA"]]$counts)
  # Name each list element with the corresponding file name
  names(data_list) <- files
  
  return(data_list)
}
test_sce <- function(data.sca){
  clust_1 <- data.sca@int_metadata$c[1]
  clust_2 <- data.sca@int_metadata$c[2]
  #filter data to only use cells types of interest
  data.sca_filt <- data.sca[,colData(data.sca)$ident == toString(clust_1) | colData(data.sca)$ident == toString(clust_2) ]
  zlm.output <- zlm(~ident, data.sca_filt)
  #zlm.lr_clust <- lrTest(zlm.output, "batch")
  #names <- rownames(zlm.lr_clust)
  #values <- zlm.lr_clust[,3,3]
  print(data.sca@int_metadata$name)
  #res <- data.frame(values,row.names = names)
  return(zlm.output)
}
add_id_cl_col <- function(data,cl){
  colData(data)["ident_select"] <- with(colData(data), ifelse(ident==toString(cl),1,0))
  # colData(data)["ident_select"] <- data.frame(colData(data)) %>% 
  #   mutate("is_ident" = as.factor(if_else(ident==toString(cl),1,0))) %>% 
  #   select("is_ident")
  
  return(data)
}

get_lfc <- function(data){
  fit.zlm <- test_sce(data)
  logfc <- getLogFC(fit.zlm)
  zlm.lr_clust <- lrTest(fit.zlm, "ident")
  values_p <- zlm.lr_clust[,3,3]
  names <- logfc$primerid
  values <- logfc$logFC
  
  if(data@int_metadata$c[1] > data@int_metadata$c[2]){
    values <- values*-1
  }
  
  
  res <- data.frame(values,values_p,names=names)
  return(res)
}


get_markers_sca <- function(data,lims = 0.6) {
  #clust1 = unlist(data@misc)[1]
  #clust2 = unlist(data@misc)[2]
  #markers <- FindMarkers(object = data, ident.1 = clust1 ,ident.2=clust2)
  markers <- data
  markers$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"<br /><br /><br />
  markers$diffexpressed[markers$values > lims & markers$values_p < (0.05/2000)] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"<br /><br /><br />
  markers$diffexpressed[markers$values < -lims & markers$values_p < (0.05/2000)] <- "DOWN"
  # Explore a bit<br /><br /><br />
  
  #head(markers[order(markers$values_p) & df$diffexpressed == 'DOWN', ])
  
  markers<- rownames_to_column(markers, var = "gene_symbol")
  markers$delabel <- ifelse(markers$gene_symbol %in% head(markers[order(markers$values_p), "gene_symbol"], 5), markers$gene_symbol, NA)
  markers <- markers %>% mutate(p_val = values_p,avg_log2FC=values) 
  return(markers)
}

data_sca <- read_data_sca(d_str,"pbmc3k-0_")
data_sca_tests <- lapply(data_sca,get_lfc)
data_sca_volc <- lapply(data_sca_tests,get_markers_sca)



volc_plot_lines <- function(meth_name,plot_name) {
  
  
  
  meth <-meth_name
  
  p <- list_rbind(data_sca_volc,names_to = "id") %>% filter(id == "orig" | id==meth) %>%
    mutate(id = ifelse(id == "orig", "orig", "new")) %>% 
    mutate(id_f = as.factor(id),lp = -log10(p_val))  %>% select(id,gene_symbol,lp,avg_log2FC,diffexpressed,delabel)%>%
    pivot_wider(names_from=id,values_from=c(lp,avg_log2FC)) %>% 
    mutate(new_diffexpressed = ifelse(!is.na(avg_log2FC_new), diffexpressed, NA)) %>% group_by(gene_symbol) %>%
    summarise(across(everything(),coalesce_by_column)) %>% mutate(diffexpressed = new_diffexpressed) %>% 
    filter(!is.na(diffexpressed)) %>% 
    #mutate(diffexpressed=data_sca_volc$bbknn$diffexpressed,delabel= data_sca_volc$bbknn$delabel) %>% 
    ggplot(aes(x = avg_log2FC_new, y = lp_new, col = diffexpressed, label = delabel,group=gene_symbol)) + 
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05/2000), col = "gray", linetype = 'dashed') + 
    geom_link(data=. %>% filter(diffexpressed != "NO"),aes(x = avg_log2FC_orig, y = lp_orig, xend = avg_log2FC_new, yend = lp_new,alpha = after_stat(index)),arrow=arrow(angle=15,length=unit(0.1,"inches")))+#,show.legend = FALSE) +
    geom_point(data= . %>% filter(diffexpressed != "NO"),aes(x = avg_log2FC_orig, y = lp_orig),alpha = 0.2) +
    geom_point(aes(x = avg_log2FC_new, y = lp_new))+
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                       labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    #scale_alpha_manual(values=c(0.1,0.1,1,1) ,breaks=c(0.5), na.value = NA) +
    scale_alpha(name = "Data type",breaks=c(0,1),labels=c("Original","Batch corrected")) +
    coord_cartesian(ylim = c(0, 300), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Genes', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    scale_x_continuous(breaks = seq(-13, 13, 2)) + # to customise the breaks in the x axis
    ggtitle(paste("",plot_name,sep=""))  # Plot title 
  #geom_text_repel(max.overlaps = Inf,box.padding = 0.5,show.legend = FALSE) # To show all labels
  p <- p + theme(plot.title = element_text(size=20, face="bold", 
                                           margin = margin(10, 0, 10, 0)),
                 panel.background = element_rect(fill = "white", 
                                                 colour = NA), panel.border = element_rect(fill = NA, 
                                                                                           colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
                 panel.grid.minor = element_line(linewidth = rel(0.5)), 
                 strip.background = element_rect(fill = "grey85", colour = "grey20"), 
                 legend.key = element_rect(fill = "white",colour = NA))
  return(p)
}

pl2 <- volc_plot_lines("bbknn","BBKNN")
pl3 <- volc_plot_lines("combat","Combat")
pl4 <- volc_plot_lines("combatseq","CombatSeq")
pl5 <- volc_plot_lines("harmony","Harmony")
pl6 <- volc_plot_lines("liger","LIGER")
pl7 <- volc_plot_lines("mnn","MNN")
pl8 <- volc_plot_lines("scvi","SCVI")
pl9 <- volc_plot_lines("seurat","Seurat")
pl10 <- volc_plot_lines("seurat_v2","Seurat - Alternative")




custom_tags = list(c("BBKNN","Combat","CombatSeq","Harmony","LIGER","MNN","SCVI","Seurat","Seurat Alternative"))

plot_all <- pl2 + pl3 + pl4 +pl5 + pl6 + pl7 + pl8 + pl9  +  plot_layout(ncol = 4,guides = "collect") &theme(legend.position = 'bottom',legend.title = element_text(size = 18), 
                                                                                                             legend.text  = element_text(size = 15),
                                                                                                             legend.key.size = unit(1, "lines"))
#plot_all <- pl1 + pl2 + pl3 + pl4 +pl5 + pl6 + pl7 + pl8 + pl9 + pl10 + plot_annotation(tag_levels = custom_tags) + plot_layout(ncol = 5,guides = "collect") & theme(plot.title = element_blank())
ggsave("plots/volcano_mast.png", plot = plot_all, width = 16, height = 8, units = "in")



#SIDE DENSITY VOLCANO PLOTS

plot_side_density <- function(data_name,plot_name) {
  #plot_name="LIGER"
  p1 <- list_rbind(data_sca_volc,names_to = "id") %>% filter(id == "orig" | id==data_name) %>%
    mutate(id = ifelse(id == "orig", "orig", "new")) %>% 
    mutate(id_f = as.factor(id),lp = -log10(p_val))  %>% select(id,gene_symbol,lp,avg_log2FC,diffexpressed,delabel,p_val)%>%
    filter(p_val>= 0.05/2000) %>% 
    ggplot( aes(x = avg_log2FC, y = lp, label = delabel,color=id)) +
    geom_vline(xintercept = c(-0.3, 0.3), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05/2000), col = "gray", linetype = 'dashed') + 
    geom_point(size = 0.5,alpha=0.5) + 
    
    geom_xsidedensity(aes(y = after_stat(density)),show.legend = FALSE,alpha=0.3) +
    geom_ysidedensity(aes(x = after_stat(density)),show.legend = FALSE,alpha=0.3) +
    
    scale_color_manual(values = c("#00AFBB", "#bb0c00"), # to set the colours of our variable  
                       labels = c("Corrected",  "Original")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = c(0, 5), xlim = c(-.5, 0.5)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Genes', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),subtitle=paste("","Low signal",sep="")) + 
    #scale_x_continuous(breaks = seq(-1, 11, 2)) + # to customise the breaks in the x axis
    #ggtitle() + # Plot title 
    geom_text_repel(max.overlaps = Inf,show.legend = FALSE) # To show all labels
  p1 <- p1 + theme(legend.position = "bottom",plot.subtitle = element_text(size=8,hjust = 0.5))+ theme_ggside_void()
  #p1 <-ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE)
  p2 <- list_rbind(data_sca_volc,names_to = "id") %>% filter(id == "orig" | id==data_name) %>%
    mutate(id = ifelse(id == "orig", "orig", "new")) %>% 
    mutate(id_f = as.factor(id),lp = -log10(p_val))  %>% select(id,gene_symbol,lp,avg_log2FC,diffexpressed,delabel,p_val)%>%
    filter(diffexpressed=="NO",p_val< 0.05/2000) %>% 
    ggplot( aes(x = avg_log2FC, y = lp, label = delabel,color=id)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05/2000), col = "gray", linetype = 'dashed') + 
    geom_xsidedensity(aes(y = after_stat(density)),show.legend = FALSE,alpha=0.3) +
    geom_ysidedensity(aes(x = after_stat(density)),show.legend = FALSE,alpha=0.3) +
    #geom_line() +
    #ggforce::geom_link2(aes(color = lp), size = 5, n = 500, lineend = "round")+
    geom_point(size = 0.5,alpha=0.5) + 
    scale_color_manual(values = c("#00AFBB", "#bb0c00"), # to set the colours of our variable  
                       labels = c("Corrected",  "Original")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = c(4, 40), xlim = c(-0.6, 0.6)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Genes', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),subtitle=paste("","High signal",sep="")) + 
    #scale_x_continuous(breaks = seq(-1, 11, 2)) + # to customise the breaks in the x axis
    #ggtitle() + # Plot title 
    geom_text_repel(max.overlaps = Inf,show.legend = FALSE) # To show all labels
  p2 <- p2 + theme(legend.position = "bottom",plot.subtitle = element_text(size=8,hjust = 0.5)) + theme_ggside_void()
  
  
  #p <- (p1 + p2) 
  
  #p<- p+ plot_layout(guides = "collect") + plot_annotation(plot_name) & theme(legend.position = 'bottom') 
  return(list("plot1" = p1, "plot2" = p2))
}

pl2 <- plot_side_density("bbknn","BBKNN")
pl3 <- plot_side_density("combat","Combat")
pl4 <- plot_side_density("combatseq","Combat-seq")
pl5 <- plot_side_density("harmony","Harmony")
pl6 <- plot_side_density("liger","LIGER")
pl7 <- plot_side_density("mnn","MNN")
pl8 <- plot_side_density("scvi","SCVI")
pl9 <- plot_side_density("seurat","Seurat")
pl10 <- plot_side_density("seurat_v2","Seurat - Alternative")


plot_all1 <- (pl2$plot1 | pl2$plot2)/  (pl3$plot1 | pl3$plot2) +plot_layout(ncol = 1,guides = "collect") +
  plot_annotation(tag_level=list(c("BBKNN","","Combat","")))&theme(legend.position = 'bottom',plot.tag.position = c(0, 1),
                                                                   plot.tag = element_text(size = 11, hjust = 0, vjust = 0))
ggsave("plots/side_density_1.png", plot = plot_all1, width = 7, height = 6, units = "in")
#plot_all2 <- pl5 + pl6 + pl7 + pl8 + pl9 + pl10 +  plot_layout(ncol = 1,guides = "collect") &theme(legend.position = 'bottom',legend.title = element_text(size = 15),
plot_all2 <- (pl4$plot1 | pl4$plot2)/  (pl5$plot1 | pl5$plot2) +plot_layout(ncol = 1,guides = "collect") +
  plot_annotation(tag_level=list(c("Combat-seq","","Harmony","")))&theme(legend.position = 'bottom',plot.tag.position = c(0, 1),
                                                                         plot.tag = element_text(size = 11, hjust = 0, vjust = 0))
ggsave("plots/side_density_2.png", plot = plot_all2, width = 7, height = 6, units = "in")

plot_all3 <- (pl6$plot1 | pl6$plot2)/  (pl7$plot1 | pl7$plot2) +plot_layout(ncol = 1,guides = "collect") +
  plot_annotation(tag_level=list(c("LIGER","","MNN","")))&theme(legend.position = 'bottom',plot.tag.position = c(0, 1),
                                                                plot.tag = element_text(size = 11, hjust = 0, vjust = 0))
ggsave("plots/side_density_3.png", plot = plot_all3, width = 7, height = 6, units = "in")

plot_all4 <- (pl8$plot1 | pl8$plot2)/  (pl9$plot1 | pl9$plot2) +plot_layout(ncol = 1,guides = "collect") +
  plot_annotation(tag_level=list(c("SCVI","","Seurat","")))&theme(legend.position = 'bottom',plot.tag.position = c(0, 1),
                                                                  plot.tag = element_text(size = 11, hjust = 0, vjust = 0))
ggsave("plots/side_density_4.png", plot = plot_all4, width = 7, height = 6, units = "in")


#Confusion matrix plot

theme_set(theme_bw())
library(reshape2)

rearrange_matrix <- function(mat,rat_mat=NULL) {
  n <- nrow(mat)
  k <- ncol(mat)
  
  # For each column, find the row index of the maximum element
  max_positions <- apply(mat, 2, which.max)
  
  # Create a data frame that keeps track of original column index and max row
  col_info <- data.frame(
    col_index = seq_len(k),
    max_row = max_positions
  )
  
  # Sort columns by their max_row
  col_info <- col_info[order(col_info$max_row), ]
  
  # We'll assign new column positions based on these sorted columns.
  # For each unique max_row, we try to place the column at column position = max_row.
  # If that spot is taken or out of range, we place it in the next free spot to the right.
  
  new_positions <- integer(k)  # to store the new positions of each column
  used_positions <- rep(FALSE, k)  # track which positions are occupied
  
  for (r in unique(col_info$max_row)) {
    # Get all columns that have max_row = r
    subset_rows <- col_info[col_info$max_row == r, ]
    
    # Start placing at column position = r, if possible
    # If r > k (more rows than columns), we start at the last column and move forward (though no "forward" beyond k)
    start_pos <- min(r, k)
    
    # If start_pos is already taken, move to the next free slot
    while (start_pos <= k && used_positions[start_pos]) {
      start_pos <- start_pos + 1
    }
    
    # Assign positions to these columns
    for (j in seq_len(nrow(subset_rows))) {
      pos <- start_pos + (j - 1)
      # Move forward if necessary to find a free spot
      while (pos <= k && used_positions[pos]) {
        pos <- pos + 1
      }
      
      if (pos > k) {
        stop("Not enough column positions to place all columns according to the rule.")
      }
      
      used_positions[pos] <- TRUE
      new_positions[subset_rows$col_index[j]] <- pos
    }
  }
  
  # Reorder the columns of the matrix according to new_positions
  if(is.null(rat_mat)){
    mat_reordered <- mat[, order(new_positions), drop = FALSE]
  }
  else{
    rat_reordered <- rat_mat[, order(new_positions), drop = FALSE]
    mat_reordered <- mat[, order(new_positions), drop = FALSE]
    return(list(mat_reordered,rat_reordered))
  }
  
  
  return(mat_reordered)
}

d_str = c("bbknn","combat","combatseq","harmony","liger","mnn","scvi","seurat")
alt_str = c("liger_v2","seurat_v2")
#d_str_neuro = d_str = c("bbknn","combat","combatseq","harmony","liger","mnn","scvi","seurat")
#rat_str = c("bbknn","combat","combatseq","harmony","liger","mnn","scvi","seurat")
read_files <- function(data,neuro=FALSE,ratio=FALSE) {
  if(neuro){
    dat = read_csv(file = paste("data/neuro_",data,"_clustcon.csv",sep = ""),col_names=F)
  }
  else if(ratio){
    #print(data)
    dat = read_csv(file = paste("data/",data,"_clustrat.csv",sep = ""),col_names=F)
  }
  else{
    print(data)
    dat = read_csv(file = paste("data/",data,"_clustcon.csv",sep = ""),col_names=F)
  }
  return(dat)
}
data_list = lapply(d_str, read_files)
data_list_alt <- lapply(alt_str, read_files)
data_list_neuro = lapply(d_str, read_files,neuro=TRUE)
data_list_rat = lapply(d_str, read_files,ratio=TRUE)

data_list_neuro <- lapply(data_list_neuro,rearrange_matrix)
data_list_alt <- lapply(data_list_alt,rearrange_matrix)

dd <- map2(data_list,data_list_rat,rearrange_matrix) %>% transpose()

data_list <- dd[[1]]
data_list_rat <- dd[[2]]
names(data_list ) <- d_str
names(data_list_neuro) <- d_str
names(data_list_rat) <- d_str
names(data_list_alt) <- alt_str

create_plot <- function(plotdata, name,neuro=FALSE) {
  data <- plotdata
  colnames(data) <- paste0("", 1:ncol(data), "")
  if(neuro){
    rownames(data) <- paste0("", 1:nrow(data), "")
  }
  else{
    rownames(data) <-c("Memory CD4","Naive - CD4","CD14+ - Mono", "B" ,"CD8+ T","FCGR3A+ Mono","NK","Dendritic","Megakaryocytes")
    
  }
  #data <- melt(data)
  
  #data <- data %>% 
  #data <- melt(data) %>% mutate(Var1 = as.factor(Var1),Var2 = as.factor(Var2)) 
  #return(data)
  p <- data %>% rownames_to_column(var="Var1") %>% pivot_longer(!Var1,names_to = "Var2",values_to = "value") %>% mutate(value = ifelse(value==0,NA,value)) %>% mutate(Var1 = factor(Var1,levels=rownames(data)),Var2 = factor(Var2,levels=colnames(data))) %>%  ggplot( aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() 
  if(neuro){
    p <- p + geom_text(aes(label = round(value, 1)),size=1.1) 
  }
  else{
    p <- p + geom_text(aes(label = round(value, 1)),size=1.8) 
  }
  
  p <- p + scale_x_discrete(position = "top") +
    scale_y_discrete(limits=rev) 
  #scale_fill_gradient(low="white", high="blue") +
  
  if(neuro){
    p <- p + scale_fill_distiller(palette = "PuRd", direction=+1, limits=c(0,1990),na.value="white") 
  }
  else{
    p <- p + scale_fill_distiller(palette = "PuRd", direction=+1, limits=c(0,600),na.value="white") 
  }
  p <- p + ggtitle(name) + xlab("After") + ylab("Before")+ 
    #p <- p + xlab("After") + ylab("Before")+ 
    guides(fill = guide_colourbar(title = "Cells")) + theme(plot.title = element_text(hjust = 0.5))
  return(p)
}
test_arrays <- function(arr,neuro=FALSE){
  k = 1
  ans =TRUE
  if(neuro){
    cell_sum = 19525
  }
  else{
    cell_sum = 2638
  }
  for (i in arr) {
    #rownames(i) <- c("Memory CD4","Naive - CD4","CD14+ - Mono", "B" ,"CD8+ T","FCGR3A+ Mono","NK","Dendritic","Megakaryocytes")
    if(sum(i)!= cell_sum){
      ans = FALSE
      print(k)
      print(sum(i))
      print(i)
    }
    else{
      ans = TRUE
    }
    if(!ans){
      return(FALSE)
    }
    k=k+1
  }
  return(TRUE)
}

return_data <- function(data){
  
}

create_ratio_plot <- function(plotdata,helpdata, name,neuro=FALSE) {
  data <- plotdata
  data[data == "NaN"] <- NA
  #data <- data %>% mutate_all(~ifelse(is.nan(.), NA, .))
  colnames(data) <- paste0("", 1:ncol(data), "")
  if(neuro){
    rownames(data) <- paste0("", 1:nrow(data), "")
  }
  else{
    rownames(data) <-c("Memory CD4","Naive - CD4","CD14+ - Mono", "B" ,"CD8+ T","FCGR3A+ Mono","NK","Dendritic","Megakaryocytes")
    
  }
  #data <- melt(data) %>% mutate(value3 = 2*abs(0.5-value))
  #data$value2 <- ifelse(melt(helpdata)$value < 10,0,1)
  #value2 = ifelse(pivot_longer(!Var1,names_to = "Var2",values_to = "value")$value<10,0,1)
  data[data == "NaN"] <- NA
  
  p <- data %>% rownames_to_column(var="Var1") %>% pivot_longer(!Var1,names_to = "Var2",values_to = "value") %>%  
    mutate(value = ifelse(value==0,NA,value)) %>% 
    mutate(Var1 = factor(Var1,levels=rownames(data)),Var2 = factor(Var2,levels=colnames(data)),value3 = 2*abs(0.5-value)) %>% 
    mutate(value2 = helpdata %>%rownames_to_column(var="Var1") %>%  
             pivot_longer(!Var1,names_to = "Var2",values_to = "value") %>% 
             mutate(dd = ifelse(value<10,0,1)) %>% pull(dd)) %>%  
    ggplot( aes(x=Var2, y=Var1, fill=value3))+ geom_tile(aes(alpha=value2)) 
  if(neuro){
    p <- p + geom_text(aes(label = round(value, 1)),size=1.1) 
  }
  else{
    p <- p + geom_text(aes(label = round(value, 1)),size=1.8) 
  }
  
  p <- p + scale_x_discrete(position = "top") + guides(alpha = "none")+
    scale_y_discrete(limits=rev)  
  #scale_fill_gradient(low="white", high="blue") +
  
  if(neuro){
    p <- p + scale_fill_distiller(palette = "YlOrRd", direction=+1, limits=c(0,1),na.value="white") 
  }
  else{
    p <- p + scale_fill_distiller(palette = "YlOrRd", direction=+1, limits=c(0,1),na.value="white") 
  }
  #p <- p + ggtitle(paste(name," - ","batch ratio",sep="")) + xlab("After") + ylab("Before")+ 
  p <- p +  xlab("After") + ylab("Before")+
    guides(fill = guide_colourbar(title = "Ratio")) + theme(plot.title = element_text(hjust = 0.5),axis.text.y=element_blank(),
                                                            panel.grid.major = element_blank(), 
                                                            panel.grid.minor =element_blank())
  return(p)
  
}

p1 <- create_plot(data_list$bbknn , "BBKNN")



p2 <- create_plot(data_list$combat , "Combat")
p3 <- create_plot(data_list$combatseq , "Combat-seq")


p4 <- create_plot(data_list$harmony , "Harmony")

p_alt_1<- create_plot(data_list$liger %>% relocate(X1),"LIGER")



p5 <- create_plot(data_list_alt$liger_v2 %>% relocate(X2,.after=X4),"LIGER alternative")

p6 <- create_plot(data_list$mnn, "MNN")


p7 <- create_plot(data_list$scvi, "SCVI")


p8 <- create_plot(data_list$seurat, "Seurat")


p_alt_2 <- create_plot(data_list_alt$seurat_v2, "Seurat alternative")

p1_n <- create_plot(data_list_neuro$bbknn, "BBKNN",TRUE)



p2_n <- create_plot(data_list_neuro$combat, "Combat",TRUE)
p3_n <- create_plot(data_list_neuro$combatseq, "Combat-seq",TRUE)

p4_n <- create_plot(data_list_neuro$harmony, "Harmony",TRUE)

p5_n <- create_plot(data_list_neuro$liger %>% relocate(X6,.after=X8), "LIGER",TRUE)



p6_n <- create_plot(data_list_neuro$mnn, "MNN",TRUE)


p7_n <- create_plot(data_list_neuro$scvi, "SCVI",TRUE)


p8_n <- create_plot(data_list_neuro$seurat %>% relocate(X2), "Seurat",TRUE)

p1_r <- create_ratio_plot(data_list_rat$bbknn,data_list$bbknn, "BBKNN")
p2_r <- create_ratio_plot(data_list_rat$combat,data_list$combat, "Combat")
p3_r <- create_ratio_plot(data_list_rat$combatseq,data_list$combatseq, "Combat-seq")
p4_r <- create_ratio_plot(data_list_rat$harmony,data_list$harmony, "Harmony")
p5_r <- create_ratio_plot(data_list_rat$liger %>%  relocate(X1),data_list$liger%>%  relocate(X1), "LIGER")
p6_r <- create_ratio_plot(data_list_rat$mnn,data_list$mnn, "MNN")
p7_r <- create_ratio_plot(data_list_rat$scvi,data_list$scvi, "SCVI")
p8_r <- create_ratio_plot(data_list_rat$seurat,data_list$seurat, "Seurat")

theme_set(theme_bw())
combined3 = ((p1 + p1_r)/ (p2 + p2_r) /(p3 + p3_r) /(p4 + p4_r)) | ((p5 + p5_r) /(p6 + p6_r) /(p7 + p7_r)/(p8 + p8_r))
tags <- c("BBKNN"," ","Combat"," ","Combat-seq","","Harmony"," ","LIGER"," ","MNN"," ","SCVI"," ","Seurat"," ")
combined3 = combined3 + plot_layout(guides="collect")+ plot_annotation(tag_levels = list(tags))&theme(plot.tag = element_text(size = 14),axis.text.y = element_text(size=8),plot.title = element_blank())

combined3
ggsave("~/vagrant/batch/batch/report/fromr/big_cc.pdf", width = 15, height = 9, units = "in")

combined = p1_n+p2_n    &  theme(legend.position = "bottom") & labs(x=NULL,y=NULL)
#combined = combined  + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
#+ plot_annotation(tag_levels = list(c("BBKNN","Combat")))
combined = combined  & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
                             axis.title=element_text(size=7,face="bold"), 
                             axis.text=element_text(size=5),
                             legend.key.height= unit(0.04, 'in'),
                             # legend.key.width= unit(0.14, 'in'),
                             plot.title=element_text(size=7))  
combined + plot_layout(guides="collect",ncol=1)
ggsave("plots/cc_table_neuro_plot1.png", width = 3, height = 5, units = "in")



combined = p3_n+p4_n    &  theme(legend.position = "bottom") & labs(x=NULL,y=NULL)
#combined = combined  + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
combined = combined   & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
                              axis.title=element_text(size=7,face="bold"), 
                              axis.text=element_text(size=5),
                              legend.key.height= unit(0.04, 'in'),
                              # legend.key.width= unit(0.14, 'in'),
                              plot.title=element_text(size=7))  
combined =combined + plot_layout(guides="collect",ncol=1)
ggsave("plots/cc_table_neuro_plot2.png", width = 3, height = 5, units = "in")



combined = p5_n+p6_n    &  theme(legend.position = "bottom") & labs(x=NULL,y=NULL)
#combined = combined  + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
combined = combined   & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
                              axis.title=element_text(size=7,face="bold"), 
                              axis.text=element_text(size=5),
                              legend.key.height= unit(0.04, 'in'),
                              # legend.key.width= unit(0.14, 'in'),
                              plot.title=element_text(size=7))  
combined = combined + plot_layout(guides="collect",ncol=1)
ggsave("plots/cc_table_neuro_plot3.png", width = 3, height = 5, units = "in")


combined = p7_n +p8_n   & theme(legend.position = "bottom") & labs(x=NULL,y=NULL)
#combined = combined  + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
combined = combined   & theme(plot.tag=element_text(size=4),legend.text=element_text(size=5),legend.title=element_text(size=5),
                              axis.title=element_text(size=7,face="bold"), 
                              axis.text=element_text(size=5),
                              legend.key.height= unit(0.04, 'in'),
                              # legend.key.width= unit(0.14, 'in'),
                              plot.title=element_text(size=7))  
combined + plot_layout(guides="collect",ncol=1)
ggsave("plots/cc_table_neuro_plot4.png", width = 3, height = 5, units = "in")