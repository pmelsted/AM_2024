library(ggplot2)
theme_set(theme_gray())
library(scales)
library(plyr)
library(tidyverse)
library(patchwork)
library(ggbreak)
library(ggfortify)
library(ggrepel)
roundUp <- function(x) 10^ceiling(log10(x))
setwd("/Users/sinant/vagrant/batch/batch")

METHODS_STR = c("BBKNN","Combat","ComBat-seq","Harmony","LIGER","MNN","SCVI","Seurat")
DATA_STR = c("pbmc3k","neuro","jejunum","pbmc4k","heart")

fix_data_plot <- function(data) {
  return (data %>%
            mutate(
              variable = case_when(
                variable == "combat" ~ "Combat",
                variable == "combatseq" ~ "ComBat-seq",
                variable == "harmony" ~ "Harmony",
                variable == "liger" ~ "LIGER",
                #variable == "ligerv2" ~ "LIGER alternative",
                variable == "mnn" ~ "MNN",
                variable == "scvi" ~ "SCVI",
                variable == "seurat" ~ "Seurat",
                variable == "seuratv2" ~ "Seurat V2",
                variable == "downsample" ~ "Downsampling",
                variable == "resample" ~ "Resampling",
                #variable == "seuratv2" ~ "Seurat alternative"
              )
            )%>% filter(variable != "Seurat V2",variable != "Resampling"))
}
create_plots <- function(files,file_str) {
  # Create the full file paths based on the directory and specified file names
  file_names = lapply(files,function(x) paste(x,file_str,"_","plot.csv",sep=""))
  file_paths <- file.path("data", file_names)
  
  # Read each specified file and store it in a list
  data_list <- lapply(file_paths, read.csv)
  
  # Name each list element with the corresponding file name
  names(data_list) <- files
  nn_data <- lapply(data_list, fix_data_plot)
  return(nn_data)
}

nn_data <- create_plots(DATA_STR,"")
nn_data_emb <- create_plots(DATA_STR,"_emb")


##CREATE THE ACTUAL GGPLOT OBJECTS

pbmc3k_plot_1 <- save_pl(nn_data$pbmc3k,'nn_rank_pbmc.png',log_scale=TRUE,ret=TRUE)
pbmc3k_plot_2 <- save_pl(nn_data_emb$pbmc3k,'nn_rank_emb_pbmc.png',log_scale=TRUE,ret=TRUE,emb=TRUE)
neuro_plot_1 <- save_pl(nn_data$neuro,'nn_rank_neuro.png',log_scale=TRUE,ret=TRUE)
neuro_plot_2 <- save_pl(nn_data_emb$neuro,'nn_rank_emb_neuro.png',log_scale=TRUE,ret=TRUE,emb=TRUE)
pbmc4k_plot1 <- save_pl(nn_data$pbmc4k, "nn_rank_pbmc4k.png",log_scale = TRUE,ret=TRUE)
pbmc4k_plot2 <- save_pl(nn_data_emb$pbmc4k, "nn_rank_emb_pbmc4k.png",log_scale = TRUE,ret=TRUE)
heart_plot_1 <- save_pl(nn_data$heart,'nn_rank_heart.png',log_scale=TRUE,ret=TRUE)
heart_plot_2 <- save_pl(nn_data_emb$heart,'nn_rank_emb_heart.png',log_scale=TRUE,ret=TRUE)
jejunum_plot_1 <- save_pl(nn_data$jejunum,'nn_rank_jejunum.png',log_scale=TRUE,ret=TRUE)
jejunum_plot_2 <- save_pl(nn_data_emb$jejunum,'nn_rank_emb_jejunum.png',log_scale=TRUE,ret=TRUE,emb=TRUE)

#heart_plot_3 <- save_pl(heart_comb,'nn_rank_comb_heart.png',log_scale=TRUE,ret=TRUE)
#read.csv("data/pbmc3k_plot.csv") %>% group_by(variable) %>% dplyr::summarise(n = n())
#simul_pbmc_plot_1 <- save_pl(nn_data$simul_pbmc,'nn_rank_simul_pbmc.png',log_scale=TRUE,ret=TRUE)
#simul_pbmc_plot_2 <- save_pl(nn_data_emb$simul_pbmc,'nn_rank_emb_simul_pbmc.png',log_scale=TRUE,ret=TRUE)
#simul_neuro_plot_1 <- save_pl(nn_data$simul_neuro,'nn_rank_simul_neuro.png',log_scale=TRUE,ret=TRUE)
#simul_neuro_plot_2 <- save_pl(nn_data_emb$simul_neuro ,'nn_rank_emb_simul_neuro.png',log_scale=TRUE,ret=TRUE)


#pbmc3k_1 <- read.csv(file = 'data/pbmc3k_plot.csv')
#pbmc3k_2 <- read.csv(file = "data/pbmc3k_emb_plot.csv")
#jejunum_1 <- read.csv(file = 'data/jejunum_plot.csv')
#jejunum_1 %>% group_by(variable) %>% summarise(n = n())



create_cc_dat <- function(){
  c_pbmc <- c(0.18716102, 0.04430804, 0.04865283, 0.23030268, 0.27210299,
              0.10315215, 0.08587508)
  c_neuro <- c(0.24135723, 0.1090653 , 0.11131882, 0.26081946, 0.32089629,
               0.15925736, 0.24847631)
  c_pbmc_simul <- c(0.16397912, 0.06454482, 0.05183122, 0.21639354, 0.26831224,
                    0.09194028, 0.07651484)
  c_neuro_simul <- c(0.28027578, 0.13032601, 0.12756427, 0.28157983, 0.33253423,
                     0.15533537, 0.23420467)
  data_names <- rep(c("PBMC3K","Neuro","PBMC3K simul","Neuro simul"),times=c(7,7,7,7))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  ratio <- c(c_pbmc,c_neuro,c_pbmc_simul,c_neuro_simul)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("bbknn", "combat","harmony","mnn","scvi","seurat","seurat v2"),4)
  return(data.frame(ratio,data_names,methods))
  
}

create_cc_dat2 <- function(){
  
  c_pbmc4k <- c(0.16679909, 0.02059228, 0.0261459 , 0.18923706, 0.17975617,
                0.07875635, 0.07742381, 0.05428041, 0.07510891)
  c_heart <- c(0.1509037 , 0.07280907, 0.08853253, 0.1776305 , 0.18727103,
               0.11135137, 0.12137727, 0.11592593, 0.09604139)
  data_names <- rep(c("PBMC4K","Heart"),times=c(9,9))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  ratio <- c(c_pbmc4k,c_heart)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("bbknn", "combat","harmony","mnn","scvi","seurat","seurat v2","downsample","resample"),2)
  return(data.frame(ratio,data_names,methods))
  
}

create_diffexp_data <- function(file_str){
  d <- read.csv(file_str,header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  c_nr_matching <- as.numeric(d[4,])
  method_l = length(c_nr_matching)
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch","Matching genes"),times=c(method_l,method_l,method_l,method_l))
  nr_genes <- c(c_clust,c_batch,c_comb,c_nr_matching)
  methods = rep(c("Original",METHODS_STR),4)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original",METHODS_STR))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}



integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

save_pl <- function(data,f_str, log_scale=FALSE,ret=FALSE ,leg=TRUE, text=TRUE,title=NULL,emb=FALSE,nofill = FALSE){
  
  p_meds <- ddply(data, .(variable), summarise, med = round(median(value),1))
  if(nofill){
    p <- ggplot(data, aes(x=variable, y=value))+ expand_limits(y=0)
  }
  else{
    p <- ggplot(data, aes(x=variable, y=value, fill=variable))+ expand_limits(y=0)
  }
  
  if(leg){
    p <- p + geom_boxplot(alpha=0.3,outlier.size = 0.1,color="grey56")
  }
  else{
    p <- p + geom_boxplot(alpha=0.3,show.legend = FALSE,outlier.size = 0.1,color="grey56")
  }
  
  # geom_text(data = p_meds, aes(x = variable, y = med, label=round(med, digits = 4)), 
  #           size = 3, vjust = -1.5) +
  
  if(text){
    p <- p +
      geom_text_repel(data = p_meds,min.segment.length = Inf, aes(x = variable, y = round(med,digits=1), label = med), 
                      #size = 3.5,bg.color="black",bg.r=.08,color="white",segment.size=0,vjust = -2.5)
                      size = 4.5,color="black",segment.size=0,vjust = -2.7,hjust=0.8)
  }
  if(!is.null(title)){
    p <- p + ggtitle(title)
  }
  # else{
  #   
  # }
  if(emb){
    p <- p + theme(axis.text.x = element_text(angle=45, vjust=.9, hjust=1))
  }
  p <- p +
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)),
          panel.background = element_rect(fill = "white", 
                                          colour = NA), panel.border = element_rect(fill = NA, 
                                                                                    colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(linewidth = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                        colour = NA)
          
    )+#ggtitle('Mean rank change of top 30 NN - Embedding')+
    xlab("Method") + ylab("Change in rank") + labs(fill = "Method") 
  if(log_scale){
    p <- p +  scale_y_continuous(trans='log10',
                                 #breaks=trans_breaks('log10', function(x) 10^x),
                                 #breaks = integer_breaks(),
                                 labels=trans_format('log10', math_format(10^.x)))
    
  } 
  
  p
  #p
  #print(paste("~/vagrant/batch/batch/report/fromr/",f_str,sep=""))
  #png(paste("~/vagrant/batch/batch/report/fromr/",f_str,sep=""))
  if(ret){
    return(p)
  }
  else{
    ggsave(paste("~/vagrant/batch/batch/report/fromr/",f_str,sep=""),width= 5.7,height = 3)
  }
  
}



##FIX NAMES IN DATA

pbmc_alt_2 <- read.csv("data/pbmc3k_emb_plot.csv") %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "liger" ~ "LIGER",
      variable == "ligerv2" ~ "LIGER alternative",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat alternative"
    )
    
  ) %>% filter(variable == "LIGER alternative"|variable =="LIGER"|variable =="Seurat alternative"|variable =="Seurat")









#PLOTTING FUNCTIONS FOR DIFFEXP PLOTS


run_plots <- function(){
  diffexp_dat <- create_diffexp_data("data/diffexp_pbmc_full.csv")
  diffexp_dat_orig <- create_diffexp_data("data/diffexp_pbmc_orig_x.csv")
  
  
  p0 = diffexp_dat %>% filter(data_names != "Matching genes") %>% ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    #geom_hline(yintercept=179,color="blue")+
    #ylim(0, 2000)+
    #coord_cartesian(ylim=c(0, 1000))
    theme(legend.position="bottom",plot.title = element_text(size=20, face="bold", 
                                                             margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab("Numer of DE genes") + labs(fill = "Model")
  p0 <- p0# +  expand_limits(y = c(0, 1000), )
  p1<- p0 + scale_y_break(c(250,800),ticklabels =c(800,1900),scales=c(0.2),space = 0.1)# + theme(legend.position="bottom") 
  #p1<- p0 + scale_wrap(2)# + theme(legend.position="bottom") 
  p1 <- p1 + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust=0.6,angle=90),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1,size=11),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  #p1 + scale_y_continuous(breaks = c(1900, 2000))
  p2 = diffexp_dat_orig %>% filter(data_names != "Matching genes") %>% ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    geom_text(aes(label=nr_genes), position=position_dodge(width=0.9), vjust=-0.25,size=2.2)+
    ylim(0, 250)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab(" ") + labs(fill = "Model")
  
  p2 <- p2 + theme(legend.position="bottom")
  p2 <- p2 + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1,size=11),
    plot.margin = margin(48, 5, 14, 0, "pt"),
    legend.title = element_blank()
    
  )
  
  c1 = diffexp_dat %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_dat %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_full = c1/c2
  #c1 = c(123,175,175,143.5,166,151,173)
  c1 = diffexp_dat_orig %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_dat_orig %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  #c1 = diffexp_dat_orig_neuro[2:8,1]
  #c2 = diffexp_dat_orig[2:8,1]
  nr_matching_genes_orig = c1/c2
  
  #nr_matching_genes_full = c(121/133, 163.5/164.5 , 174/177, 137/155.5,172.5/897, 87/124,
                             #176/877.5)
  #nr_matching_genes_orig =  c( 121/133,175/177,174/177, 137/155.5,161/167, 153/196.5,
                              # 172.5/179)
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(METHODS_STR,2)
  type = c(rep("Corrected counts",length(METHODS_STR)),rep("Uncorrected Counts",length(METHODS_STR)))
  full_rat = data.frame(nr_matching_genes,methods,type)
  
  
  p3 = ggplot(full_rat, aes( y=nr_matching_genes, x=as.factor(methods),fill=type)) + 
    geom_bar(alpha=0.9,stat='identity',position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", nr_matching_genes), y= nr_matching_genes),  vjust = 1.5,
              position =position_dodge(width = 0.9))+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Method") + ylab("Percentage of matching genes") + labs(fill = "Method") + scale_fill_brewer(palette="Dark2")
  
  p3 <- p3 + theme(
    axis.text.x = element_text(size = 11,angle=30),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    #axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  return(list(p1,p2,p3))
}
p_list <-run_plots()
p1 <- p_list[1][[1]]
p2 <- p_list[2][[1]]
p3 <- p_list[3][[1]]

run_plots_neuro<- function(){
  diffexp_dat_neuro <- create_diffexp_data("data/diffexp_neuro_full.csv")
  diffexp_dat_orig_neuro <- create_diffexp_data("data/diffexp_neuro_orig_x.csv")
  
  p0_neuro = diffexp_dat_neuro %>%filter(data_names != "Matching genes") %>%  ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    #geom_hline(yintercept=831,color="blue")+
    #ylim(0, 1700)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab("Numer of DE genes") + labs(fill = "Model")
  p1_neuro<- p0_neuro + scale_y_break(c(1000,1200),ticklabels =c(1200,1500),scales=c(0.2),space = 0.1)# + theme(legend.position="bottom") 
  #p1<- p0 + scale_wrap(2)# + theme(legend.position="bottom") 
  p1_neuro <- p1_neuro + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust=0.6,angle=90),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1,size=11),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  p2_neuro = diffexp_dat_orig_neuro %>% filter(data_names != "Matching genes") %>%ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    ylim(0, 1000)+
    geom_text(aes(label=nr_genes), position=position_dodge(width=0.9), vjust=-0.25,size=2.2)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab(" ") + labs(fill = "Model")
  
  p2_neuro <- p2_neuro + theme(legend.position="bottom")
  p2_neuro <- p2_neuro + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1,size=11),
    plot.margin = margin(48, 5, 14, 0, "pt"),
    legend.title = element_blank()
    
  )
  

  
  #c1 = c(386,812.5,817,365.5,753.5,348.5,800.5)
  c1 = diffexp_dat_neuro %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_dat_neuro %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_full = c1/c2
  #c1 = c(386,819,817,365.5,798,378,820.5)
  c1 = diffexp_dat_orig_neuro %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_dat_orig_neuro %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_orig = c1/c2
  
  
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(METHODS_STR,2)
  type = c(rep("Corrected counts",length(METHODS_STR)),rep("Uncorrected Counts",length(METHODS_STR)))
  full_rat = data.frame(nr_matching_genes,methods,type)
  #methods = c("bbknn", "combat","harmony","mnn","scvi","seurat")
  #full_rat = data.frame(nr_matching_genes_neuro,methods)
  p3_neuro = ggplot(full_rat, aes( y=nr_matching_genes, x=as.factor(methods),fill=type)) + 
    geom_bar(alpha=0.9,stat='identity',position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", nr_matching_genes), y= nr_matching_genes),  vjust = 3,
              position =position_dodge(width = 0.9))+
    
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Method") + ylab("Percentage of matching genes") + labs(fill = "Method") + scale_fill_brewer(palette="Dark2")
  
  p3_neuro <- p3_neuro + theme(
    axis.text.x = element_text(size = 11,angle=30),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    #axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  return(list(p1_neuro,p2_neuro,p3_neuro))
}
p_list <-run_plots_neuro()
p1_neuro <- p_list[1][[1]]
p2_neuro <- p_list[2][[1]]
p3_neuro <- p_list[3][[1]]



run_plots_simul_pbmc<- function(){
  diffexp_simul_pbmc <- create_diffexp_data("data/diffexp_simul_pbmc_full.csv")
  diffexp_simul_pbmc_orig <- create_diffexp_data("data/diffexp_simul_pbmc_orig_x.csv")
  p0_simul_pbmc = diffexp_simul_pbmc %>% filter(data_names != "Matching genes") %>% ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_hline(yintercept=161,color="blue")+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    #ylim(0, 2000)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab("Numer of genes") + labs(fill = "Model")
  
  #p1_simul_pbmc <- p0_simul_pbmc# +  expand_limits(y = c(0, 1000), )
  p1_simul_pbmc<- p0_simul_pbmc + scale_y_break(c(250,800),ticklabels =c(800,1900),scales=c(0.2),space = 0.1)# +
  p1_simul_pbmc <- p1_simul_pbmc + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust=0.6,angle=90),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  
  p2_simul_pbmc = diffexp_simul_pbmc_orig %>% filter(data_names != "Matching genes")  %>% ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    geom_text(aes(label=nr_genes), position=position_dodge(width=0.9), vjust=-0.25,size=2.2)+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    ylim(0, 250)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab(" ") + labs(fill = "Model")
  
  p2_simul_pbmc <- p2_simul_pbmc + theme(legend.position="bottom")
  p2_simul_pbmc <- p2_simul_pbmc + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    plot.margin = margin(39, 5, 14, 0, "pt"),
    legend.title = element_blank()
    
  )
  #c1 = c(114,144,156,122,156,124,158)
  c1 = diffexp_simul_pbmc %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_simul_pbmc %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_full = c1/c2
  #nr_matching_genes_full
  #c1 = c(114,156,156,122,149,142,155)
  c1 = diffexp_simul_pbmc_orig %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_simul_pbmc_orig %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_orig = c1/c2
  #nr_matching_genes_orig
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(METHODS_STR,2)
  type = c(rep("Corrected counts",length(METHODS_STR)),rep("Uncorrected Counts",length(METHODS_STR)))
  full_rat = data.frame(nr_matching_genes,methods,type)
  #methods = c("bbknn", "combat","harmony","mnn","scvi","seurat")
  #full_rat = data.frame(nr_matching_genes_neuro,methods)
  p3_simul_pbmc = ggplot(full_rat, aes( y=nr_matching_genes, x=as.factor(methods),fill=type)) + 
    geom_bar(alpha=0.9,stat='identity',position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", nr_matching_genes), y= nr_matching_genes),  vjust = 1.5,
              position =position_dodge(width = 0.9))+
    
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Method") + ylab("Percentage of matching genes") + labs(fill = "Method")+ scale_fill_brewer(palette="Dark2")
  
  p3_simul_pbmc <- p3_simul_pbmc + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    #axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  
  return(list(p1_simul_pbmc,p2_simul_pbmc,p3_simul_pbmc))
}
p_list <-run_plots_simul_pbmc()
p1_simul_pbmc <- p_list[1][[1]]
p2_simul_pbmc <- p_list[2][[1]]
p3_simul_pbmc <- p_list[3][[1]]

run_plots_simul_neuro<- function(){
  diffexp_simul_neuro <- create_diffexp_data("data/diffexp_simul_neuro_full.csv")
  diffexp_simul_neuro_orig <- create_diffexp_data("data/diffexp_simul_neuro_orig_x.csv")
  p1_simul_neuro =diffexp_simul_neuro %>% filter(data_names != "Matching genes") %>%  ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    #geom_hline(yintercept=828,color="blue")+
    ylim(0, 1650)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab("Numer of genes") + labs(fill = "Model")
  
  #p1_simul_neuro<- p1_simul_neuro + scale_y_break(c(1100,1200),ticklabels =c(1200,1500),scales=c(0.2),space = 0.1)# +
  p1_simul_neuro <- p1_simul_neuro + theme(
    axis.title.x = element_blank(),
    #axis.title.y = element_text(hjust=0.6,angle=90),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  p2_simul_neuro = diffexp_simul_neuro_orig %>% filter(data_names != "Matching genes") %>% ggplot(aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    geom_text(aes(label=nr_genes), position=position_dodge(width=0.9), vjust=-0.25,size=2.2)+
    ylim(0, 1650)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab(" ") + labs(fill = "Model")
  
  p2_simul_neuro <- p2_simul_neuro + theme(legend.position="bottom")
  p2_simul_neuro <- p2_simul_neuro + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(39, 5, 14, 0, "pt"),
    legend.title = element_blank()
    
  )
  #c1 = c(387,801.5,812.5,356,750,355.5,802.5)
  c1 = diffexp_simul_neuro %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_simul_neuro %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_full = c1/c2
  #nr_matching_genes_full
  #c1 = c(387,810.5,812.5,356.5,796.5,381.5,811.5)
  c1 = diffexp_simul_neuro_orig %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_simul_neuro_orig %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_orig = c1/c2
  #nr_matching_genes_orig
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(METHODS_STR,2)
  type = c(rep("Corrected counts",length(METHODS_STR)),rep("Uncorrected Counts",length(METHODS_STR)))
  full_rat = data.frame(nr_matching_genes,methods,type)
  
  
  
  p3_simul_neuro= ggplot(full_rat, aes( y=nr_matching_genes, x=as.factor(methods),fill=type)) + 
    geom_bar(alpha=0.9,stat='identity',position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", nr_matching_genes), y= nr_matching_genes),  vjust = 3,
              position =position_dodge(width = 0.9))+
    
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Method") + ylab("Percentage of matching genes") + labs(fill = "Method")+ scale_fill_brewer(palette="Dark2")
  
  
  p3_simul_neuro <- p3_simul_neuro + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    #axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  return(list(p1_simul_neuro,p2_simul_neuro,p3_simul_neuro))
}
p_list <-run_plots_simul_neuro()
p1_simul_neuro <- p_list[1][[1]]
p2_simul_neuro <- p_list[2][[1]]
p3_simul_neuro <- p_list[3][[1]]


run_plots_jejunum <- function(){
  diffexp_dat_jejunum <- create_diffexp_data("data/diffexp_jejunum_full.csv")
  diffexp_dat_orig_jejunum <- create_diffexp_data("data/diffexp_jejunum_orig_x.csv")
  
  
  p0 = diffexp_dat_jejunum %>% filter(data_names != "Matching genes") %>% ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    #geom_hline(yintercept=179,color="blue")+
    #ylim(0, 2000)+
    #coord_cartesian(ylim=c(0, 1000))
    theme(legend.position="bottom",plot.title = element_text(size=20, face="bold", 
                                                             margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab("Numer of DE genes") + labs(fill = "Model")
  p0 <- p0# +  expand_limits(y = c(0, 1000), )
  p1<- p0 + scale_y_break(c(300,700),ticklabels =c(700,1900),scales=c(0.2),space = 0.1)# + theme(legend.position="bottom") 
  #p1<- p0 + scale_wrap(2)# + theme(legend.position="bottom") 
  p1 <- p1 + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust=0.6,angle=90),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1,size=11),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  #p1 + scale_y_continuous(breaks = c(1900, 2000))
  p2 = diffexp_dat_orig_jejunum %>% filter(data_names != "Matching genes") %>% ggplot( aes( y=nr_genes, x=methods,fill=data_names)) + 
    #geom_bar(alpha=0.7,stat='identity') +
    #geom_line() +
    #geom_point() +
    geom_bar(stat="identity",position=position_dodge())+
    #geom_text(aes(label = sprintf("%.1f", mean_cc_change), y= mean_cc_change),  vjust = 3)+
    geom_text(aes(label=nr_genes), position=position_dodge(width=0.9), vjust=-0.25,size=2.2)+
    ylim(0, 300)+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#,axis.text.x=element_text(angle=45,vjust=1,hjust=1))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Data") + ylab(" ") + labs(fill = "Model")
  
  p2 <- p2 + theme(legend.position="bottom")
  p2 <- p2 + theme(
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1,size=11),
    plot.margin = margin(t =34, r=5, b=14, l=0, "pt"),
    legend.title = element_blank()
    
  )
  
  c1 = diffexp_dat_jejunum %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_dat_jejunum %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  nr_matching_genes_full = c1/c2
  #c1 = c(123,175,175,143.5,166,151,173)
  c1 = diffexp_dat_orig_jejunum %>% filter(data_names == "Matching genes") %>% select(nr_genes) %>% slice(-1) %>% pull()
  c2 = diffexp_dat_orig_jejunum %>% filter(data_names == "Cluster only") %>% select(nr_genes) %>% slice(-1) %>% pull()
  #c1 = diffexp_dat_orig_neuro[2:8,1]
  #c2 = diffexp_dat_orig[2:8,1]
  nr_matching_genes_orig = c1/c2
  
  #nr_matching_genes_full = c(121/133, 163.5/164.5 , 174/177, 137/155.5,172.5/897, 87/124,
  #176/877.5)
  #nr_matching_genes_orig =  c( 121/133,175/177,174/177, 137/155.5,161/167, 153/196.5,
  # 172.5/179)
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(METHODS_STR,2)
  type = c(rep("Corrected counts",length(METHODS_STR)),rep("Uncorrected Counts",length(METHODS_STR)))
  full_rat = data.frame(nr_matching_genes,methods,type)
  
  
  p3 = ggplot(full_rat, aes( y=nr_matching_genes, x=as.factor(methods),fill=type)) + 
    geom_bar(alpha=0.9,stat='identity',position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", nr_matching_genes), y= nr_matching_genes),  vjust = 1.5,
              position =position_dodge(width = 0.9))+
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Method") + ylab("Percentage of matching genes") + labs(fill = "Method") + scale_fill_brewer(palette="Dark2")
  
  p3 <- p3 + theme(
    axis.text.x = element_text(size = 11,angle=30),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA), complete = TRUE,legend.position="bottom",
    #axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    #plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  return(list(p1,p2,p3))
}
p_list_jejunum <-run_plots_jejunum()
p1_jejunum <- p_list_jejunum[1][[1]]
p2_jejunum <- p_list_jejunum[2][[1]]
p3_jejunum <- p_list_jejunum[3][[1]]


#PATCHWORK OBJECTS

#DE plot

combined_de = ((p1 + p2)/p3) | ((p1_neuro + p2_neuro)/p3_neuro)
#theme_get()$plot.margin 
combined_de = combined_de + plot_layout(guides="collect") + plot_annotation(tag_levels = "A") & 
  theme(plot.tag=element_text(size=9),legend.position = "bottom",plot.margin=margin(r=7.5,l=5.5,b=5.5,t=5.5))
combined_de
ggsave("plots/de_combined.pdf", width = 15.5, height = 9, units = "in")

#,axis.text.x = element_text(size = 14) 

#DE - simulations plot pmbc3k

combined_simul_pbmc = p1_simul_pbmc + p2_simul_pbmc + wrap_elements(full = p3_simul_pbmc)  #& theme(legend.position = "bottom") & labs(x=NULL)
#combined = (p1 | p2)/p3 & theme(legend.position = "bottom") & labs(x=NULL)
combined_simul_pbmc <- combined_simul_pbmc + plot_annotation(tag_levels = "A")  & theme(plot.tag=element_text(size=9),legend.position = "bottom" )
layout <- "
AABB
CCCC
"
combined_simul_pbmc + plot_layout(design=layout)
ggsave("plots/de_simul_pbmc.pdf", width = 8, height = 7, units = "in")



#DE - simulations plot neuro


combined_simul_neuro = p1_simul_neuro + p2_simul_neuro+ wrap_elements(full = p3_simul_neuro) #& theme(legend.position = "bottom") & labs(x=NULL)
combined_simul_neuro <- combined_simul_neuro + plot_annotation(tag_levels = "A")  & theme(plot.tag=element_text(size=9),legend.position = "bottom" )
layout <- "
AABB
CCCC
"
combined_simul_neuro + plot_layout(design=layout,heights=c(3,6))
ggsave("plots/de_simul_neuro.pdf", width = 8, height = 7, units = "in")

#DE JEJUNUM

combined_jejunum = p1_jejunum + p2_jejunum+ wrap_elements(full = p3_jejunum) #& theme(legend.position = "bottom") & labs(x=NULL)
combined_jejunum <- combined_jejunum + plot_annotation(tag_levels = "A")  & 
  theme(plot.tag=element_text(size=9),legend.position = "bottom",plot.tag.position = c(0, 1) )
layout <- "
AABB
CCCC
"
combined_jejunum+ plot_layout(design=layout)
ggsave("plots/de_jejunum.pdf", width = 8, height = 7, units = "in")

p1_jejunum
p2_jejunum
p3_jejunum


# combined_de = ((p1_simul_pbmc + p2_simul_pbmc)/p3_simul_pbmc) | ((p1_simul_neuro + p2_simul_neuro)/p3_simul_neuro)
# combined_de = combined_de + plot_layout(guides="collect") + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9),legend.position = "bottom")
# combined_de
# ggsave("plots/de_simulated_combined.png", width = 15, height = 9, units = "in")








##Neighbor plots




neuro_plot_1 = neuro_plot_1 + ylab(NULL)
combined = pbmc3k_plot_1 + neuro_plot_1 & theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) & labs(x=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(0.4,10000))
combined = combined + plot_layout(guides="collect",ncol=2) + labs(tag = "level") 
combined

ggsave("plots/nn_plot.pdf", width = 6.5, height = 4, units = "in")
##Neighbor plots - embedding
neuro_plot_2= neuro_plot_2 + ylab(NULL)
combined = pbmc3k_plot_2 + neuro_plot_2 & theme(legend.position = "none") & labs(x=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(0.4,10000))
combined = combined + plot_layout(guides="collect",ncol=2) + labs(tag = "level") 
ggsave("plots/nn_plot2.pdf", width = 6.5, height = 4, units = "in")

# Neighbor plots - other -emb

combined =  heart_plot_2 + pbmc4k_plot2& theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) & labs(x=NULL,y=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(1,4500))
combined = combined + plot_layout(guides="collect",ncol=2)
ggsave("plots/rank_other.pdf", width = 9, height = 5, units = "in")
# Neighbor plots - other -nonemb

combined =  heart_plot_1 + pbmc4k_plot1& theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) & labs(x=NULL,y=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(1,4500))
combined = combined + plot_layout(guides="collect",ncol=2)
ggsave("plots/rank_other-non-emb.pdf", width = 9, height = 5, units = "in")

alt_plot_2 <- save_pl(pbmc_alt_2,'nn_rank_pbmc.png',log_scale=TRUE,ret=TRUE,nofill=T)

ggsave("plots/rank_alt.pdf", width = 7, height = 5, units = "in")

## neighbour nn and emb - jejunum

jejunum_plot_2 = jejunum_plot_2 + ylab(NULL)
jejunum_plot_1 = jejunum_plot_1  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
combined = jejunum_plot_1 + jejunum_plot_2 & labs(x=NULL) & theme(legend.position = "none")
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(0.4,10000))
combined = combined  +plot_layout(ncol=2) + labs(tag = "level") 
combined
ggsave("plots/rank_jejunum.pdf", width = 6.5, height = 4, units = "in")

## harmony confusion run
label=c(rep("100 PCA", 26),rep("50 PCA", 26))
pca  = c(2,4,6,8,10,12,14,18,20,22,24,26,28,30,32,34,36,38,40,42,50,60,70,80,90,100) 
nn1 = c(213.0,124.5,111.0,105.5,89.0,79.5,72.0,57.5,53.0,47.5,42.5,39.0,35.0,32.5,29.5,27.0,25.0,23.5,21.5,20.5,16.5,13.0,10.0,7.0,4.5,4)
cc1 = c(0.4584912812736922
        ,0.300606520090978
        ,0.22441243366186503
        ,0.22327520849128127
        ,0.19598180439727067
        ,0.1497346474601971
        ,0.15295678544351782
        ,0.13362395754359363
        ,0.13646702047005305
        ,0.09287338893100834
        ,0.10348749052312359
        ,0.10064442759666414
        ,0.09742228961334345
        ,0.09533737680060653
        ,0.06823351023502652
        ,0.07107657316148597
        ,0.09666413949962092
        ,0.06671721000758149
        ,0.08377558756633813
        ,0.07050796057619407
        ,0.09818043972706594
        ,0.0640636846095527
        ,0.059893858984078854
        ,0.05477634571645185
        ,0.041319181197877176,
        0.042)
nn2 = c(182.75,
        99.0,
        86.75,
        76.0,
        60.0,
        52.0,
        43.5,
        32.5,
        28.0,
        24.0,
        21.0,
        18.5,
        17.0,
        15.0,
        13.5,
        12.0,
        10.5,
        9.5,
        8.0,
        7.0,
        1.0,
        NA,
        NA,
        NA,
        NA,
        NA)
cc2= c(0.46531463229719483,
       0.28582259287338896,
       0.24374526156178922,
       0.2086808188021228,
       0.19749810462471568,
       0.1961713419257013,
       0.18233510235026534,
       0.1309704321455648,
       0.12869598180439729,
       0.10538286580742987,
       0.08188021228203185,
       0.12623199393479909,
       0.07865807429871113,
       0.08510235026535254,
       0.08927217589082638,
       0.07012888551933281,
       0.07391963608794541,
       0.05913570887035634,
       0.07354056103108415,
       0.09761182714177406,
       0.03184230477634571,
       NA,
       NA,
       NA,
       NA,
       NA)

df = data.frame(cc=c(cc1,cc2),nn=c(nn1,nn2),label,pca=pca)
p_harm = ggplot(data=df, aes(x=pca, y=cc, group=label, color=label)) +
  geom_line() + geom_point() +guides(color=guide_legend(title="Ground truth")) +
  ylab("Consensus cluster metric") + xlab("Number of PCA components") + theme(
    legend.position="bottom",
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA))

p_harm2 = ggplot(data=df, aes(x=pca, y=nn, group=label, color=label)) +
  geom_line() + geom_point() +guides(color=guide_legend(title="Ground truth")) +
  ylab("median NN rank displacement") + xlab("Number of PCA components")+ theme(
    legend.position="bottom",
    panel.background = element_rect(fill = "white", 
                                    colour = NA), panel.border = element_rect(fill = NA, 
                                                                              colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(linewidth = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                  colour = NA))
combined_harm = p_harm+ p_harm2
combined_harm = combined_harm + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9),legend.position="bottom")
combined_harm + plot_layout(guides="collect",ncol=2)
ggsave("plots/harmony_peturbation.pdf", width = 9, height = 5, units = "in")

