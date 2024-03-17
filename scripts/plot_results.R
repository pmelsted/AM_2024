theme_set(theme_gray())
library(scales)
library(plyr)
roundUp <- function(x) 10^ceiling(log10(x))




pbmc_1 <- read.csv(file = 'data/pbmc_plot.csv')
pbmc_2 <- read.csv(file = 'data/pbmc_emb_plot.csv')


neuro_1 <- read.csv(file = 'data/neuro_plot.csv')
neuro_2 <- read.csv(file = 'data/neuro_emb_plot.csv')

pbmc4k_1 <- read.csv(file="data/pbmc4k_plot.csv")
pbmc4k_2 <- read.csv(file="data/pbmc4k_emb_plot.csv")
heart_1 <-read.csv(file = 'data/heart_plot.csv')
heart_2 <- read.csv(file = 'data/heart_emb_plot.csv')


simul_pbmc_1 <- read.csv('data/simul_pbmc_plot.csv')
simul_pbmc_2 <- read.csv('data/simul_pbmc_emb_plot.csv')

simul_neuro_1 <- read.csv('data/simul_neuro_plot.csv')
simul_neuro_2 <- read.csv('data/simul_neuro_emb_plot.csv')

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

create_diffexp_dat <- function(){
  #c_clust <- c(179. , 133. , 163. , 178. , 910.5, 128. , 879.5)
  #c_batch <- c(0.0000e+00, 0.0000e+00, 1.9490e+03, 0.0000e+00, 1.8845e+03,
  #     1.0000e+00, 1.8450e+03)
  #c_comb <- c(154 ,  115 ,  1953.5 ,  153.5, 1894,  123.5, 1867)
  d <- read.csv("data/diffexp_pbmc_full.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER", "MNN","SCVI","Seurat"),3)
  #return(data.frame(nr_genes,data_names,methods))
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}

create_diffexp_dat_orig <- function(){
  #c_clust <- c(179,177.5 , 172.5 , 199 , 179)
  #c_batch <- c(0,0, 0, 0, 0)
  #c_comb <- c(178,177. ,  172.5 ,  200 ,  178.5)
  
  #c_clust <- c(179,133,177,177,167,196.5,179)
  #c_batch <- c(0,0, 0, 0, 0,0,0)
  #c_comb <- c(154,115,  153.5 ,  153.5 ,  150,172.5,156)
  d <- read.csv("data/diffexp_pbmc_orig.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER", "MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  #return(data.frame(nr_genes,data_names,methods))
  return(data)
  
}


create_diffexp_dat_neuro <- function(){
  #c_clust <- c(831 , 624 , 830 , 826 , 1036.5, 296.5 , 1021.5)
  #c_batch <- c(1, 1, 109, 1, 1317,
  #    0, 769)
  #c_comb <- c(780.5 ,  590 ,  859.5 ,  772, 1609.5,  268.5, 1318.5)
  d <- read.csv("data/diffexp_neuro_full.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}

create_diffexp_dat_orig_neuro <- function(){
  # c_clust <- c(1089,1087 , 1078.5 , 1047 , 1083)
  # c_batch <- c(1,1, 1, 1, 1)
  # c_comb <- c(1084.5,1083 ,  1074 ,  1047 ,  1077)
  #c_clust <- c(831,624 , 829 , 826 , 824,823.5,833.5)
  #c_batch <- c(1,1, 1, 1, 1,1,1)
  #c_comb <- c(780.5,590 ,  780 ,  772,  770,776,786)
  d <- read.csv("data/diffexp_neuro_orig.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER", "MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}

create_diffexp_simul_pbmc <- function(){
  #c_clust <- c(161 , 123 , 146 , 161 , 833.5, 123.5 , 875.5)
  #c_batch <- c(70.5, 67, 1925, 69,1876.5,
  #     1, 1849)
  #c_comb <- c(185 ,  152,  1931.5 ,  184, 1888.5,  121, 1869)
  d <- read.csv("data/diffexp_pbmc_full_simul.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}

create_diffexp_simul_pbmc_orig <- function(){
  #c_clust <- c(161,123 , 161 , 161 , 158,179,165)
  #c_batch <- c(70.5,67, 70, 69, 59.5,70,68)
  #c_comb <- c(185,152,  184.5 ,  184 ,  182.5,202,186)
  d <- read.csv("data/diffexp_pbmc_orig_simul.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}

create_diffexp_simul_neuro<- function(){
  # c_clust <- c(823 , 630 , 830.5 , 833.5 , 1043, 306.5 , 1096)
  # c_batch <- c(36, 19.5, 207.5, 36, 1318.5,
  #      0, 774.5)
  # c_comb <- c(779 ,  603.5,  901 ,  782, 1630.5,  277, 1422)
  d <- read.csv("data/diffexp_neuro_full_simul.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
  data <- data  %>% filter(data_names != "Batch only")
  return(data)
  
}

create_diffexp_simul_neuro_orig <- function(){
  # c_clust <- c(828,630 , 832.5 , 833.5 , 839,825,839)
  # c_batch <- c(36,19.5, 36, 36, 33.5,34.5,36)
  # c_comb <- c(779,603.5 ,  782 ,  782 ,  790.5,771.5,788)
  d <- read.csv("data/diffexp_neuro_orig_simul.csv",header=FALSE)
  c_clust <- as.numeric(d[1,])
  c_batch <- as.numeric(d[2,])
  c_comb <- as.numeric(d[3,])
  data_names <- rep(c("Cluster only","Batch only","Cluster + batch"),times=c(8,8,8))
  #meth_names <- c("PBMC3K","Neuro","PBMC4K down","Heart down","PBMC4K resamp","Heart resamp","PBMC3K simul","Neuro simul")
  nr_genes <- c(c_clust,c_batch,c_comb)
  #x = rep(c(1,2,3,4,5,6,7,8),7)
  methods = rep(c("Original","BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),3)
  data <-data.frame(nr_genes,data_names,methods)
  data$methods <- factor(data$methods,levels = c("Original","BBKNN","Combat","Harmony","LIGER","MNN","SCVI","Seurat"))
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

save_pl <- function(data,f_str, log_scale=FALSE,ret=FALSE ,leg=TRUE, text=TRUE,title=NULL,emb=FALSE){
  
  p_meds <- ddply(data, .(variable), summarise, med = round(median(value),1))
  p <- ggplot(data, aes(x=variable, y=value, fill=variable))+ expand_limits(y=0)
  
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
      geom_text_repel(data = p_meds, aes(x = variable, y = round(med,digits=1), label = med), 
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

pbmc_1 <- pbmc_1 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2"
    )
    
  ) %>% filter(variable != "Seurat V2")
neuro_1 <- neuro_1 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2"
    )
  )%>% filter(variable != "Seurat V2")

pbmc_2 <- pbmc_2 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "liger" ~ "LIGER",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2"
    )
    
  ) %>% filter(variable != "Seurat V2")

neuro_2 <- neuro_2 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "liger" ~ "LIGER",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2"
    )
  )%>% filter(variable != "Seurat V2")

heart_1 <- heart_1 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "Liger" ~ "LIGER",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2",
      variable == "downsample" ~ "Downsampling",
      variable == "resample" ~ "Resampling"
    )
  )%>% filter(variable != "Seurat V2",variable != "Resampling")

pbmc4k_1 <- pbmc4k_1 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "Liger" ~ "LIGER",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2",
      variable == "downsample" ~ "Downsampling",
      variable == "resample" ~ "Resampling"
    )
  )%>% filter(variable != "Seurat V2",variable != "Resampling")
heart_2 <- heart_2 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "liger" ~ "LIGER",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2",
      variable == "downsample" ~ "Downsampling",
      variable == "resample" ~ "Resampling"
    )
  )%>% filter(variable != "Seurat V2",variable != "Resampling")

pbmc4k_2 <- pbmc4k_2 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "liger" ~ "LIGER",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat V2",
      variable == "downsample" ~ "Downsampling",
      variable == "resample" ~ "Resampling"
    )
  )%>% filter(variable != "Seurat V2",variable != "Resampling")


pbmc_plot_1 <- save_pl(pbmc_1,'nn_rank_pbmc.png',log_scale=TRUE,ret=TRUE)
pbmc_plot_2 <- save_pl(pbmc_2,'nn_rank_emb_pbmc.png',log_scale=TRUE,ret=TRUE,emb=TRUE)
neuro_plot_1 <- save_pl(neuro_1,'nn_rank_neuro.png',log_scale=TRUE,ret=TRUE)
neuro_plot_2 <- save_pl(neuro_2,'nn_rank_emb_neuro.png',log_scale=TRUE,ret=TRUE,emb=TRUE)
pbmc4k_plot1 <- save_pl(pbmc4k_1, "nn_rank_pbmc4k.png",log_scale = TRUE,ret=TRUE)
pbmc4k_plot2 <- save_pl(pbmc4k_2, "nn_rank_emb_pbmc4k.png",log_scale = TRUE,ret=TRUE)
heart_plot_1 <- save_pl(heart_1,'nn_rank_heart.png',log_scale=TRUE,ret=TRUE)
heart_plot_2 <- save_pl(heart_2,'nn_rank_emb_heart.png',log_scale=TRUE,ret=TRUE)
#heart_plot_3 <- save_pl(heart_comb,'nn_rank_comb_heart.png',log_scale=TRUE,ret=TRUE)

simul_pbmc_plot_1 <- save_pl(simul_pbmc_1,'nn_rank_simul_pbmc.png',log_scale=TRUE,ret=TRUE)
simul_pbmc_plot_2 <- save_pl(simul_pbmc_2,'nn_rank_emb_simul_pbmc.png',log_scale=TRUE,ret=TRUE)
simul_neuro_plot_1 <- save_pl(simul_neuro_1 ,'nn_rank_simul_neuro.png',log_scale=TRUE,ret=TRUE)
simul_neuro_plot_2 <- save_pl(simul_neuro_2 ,'nn_rank_emb_simul_neuro.png',log_scale=TRUE,ret=TRUE)

pbmc_alt_2 <- pbmc_2 %>%
  mutate(
    variable = case_when(
      variable == "combat" ~ "Combat",
      variable == "harmony" ~ "Harmony",
      variable == "LIGER" ~ "LIGER",
      variable == "LIGERv2" ~ "LIGER alternative",
      variable == "mnn" ~ "MNN",
      variable == "scvi" ~ "SCVI",
      variable == "seurat" ~ "Seurat",
      variable == "seuratv2" ~ "Seurat alternative"
    )
    
  ) %>% filter(variable == "LIGER alternative"|variable =="LIGER"|variable =="Seurat alternative"|variable =="Seurat")
alt_plot_2 <- save_pl(pbmc_alt_2,'nn_rank_pbmc.png',log_scale=TRUE,ret=TRUE)
alt_plot_2


##Neighbor plots

neuro_plot_1 = neuro_plot_1 + ylab(NULL)
combined = pbmc_plot_1 + neuro_plot_1 & theme(legend.position = "bottom") & labs(x=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))
combined + plot_layout(guides="collect",ncol=2) + labs(tag = "level") 
ggsave("plots/nn_plot.png", width = 6.5, height = 4, units = "in")
##Neighbor plots - embedding
neuro_plot_2= neuro_plot_2 + ylab(NULL)
combined = pbmc_plot_2 + neuro_plot_2 & theme(legend.position = "bottom") & labs(x=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(1,10000))
combined + plot_layout(guides="collect",ncol=2) + labs(tag = "level") 
ggsave("plots/nn_plot2.png", width = 6.5, height = 4, units = "in")

# Neighbor plots - other

combined =  heart_plot_2 + pbmc4k_plot2& theme(legend.position = "bottom",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) & labs(x=NULL,y=NULL)
combined = combined + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9))&
  scale_y_continuous(trans='log10',
                     labels=trans_format('log10', math_format(10^.x)),limits=c(1,4500))
combined + plot_layout(guides="collect",ncol=2)
ggsave("plots/rank_other.png", width = 9, height = 5, units = "in")

run_plots <- function(){
  diffexp_dat <- create_diffexp_dat()
  diffexp_dat_orig <- create_diffexp_dat_orig()
  
  
  p0 = ggplot(diffexp_dat, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  #p1 + scale_y_continuous(breaks = c(1900, 2000))
  p2 = ggplot(diffexp_dat_orig, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    plot.margin = margin(48, 5, 14, 0, "pt"),
    legend.title = element_blank()
    
  )
  nr_matching_genes_full = c(121/133, 163.5/164.5 , 174/177, 137/155.5,172.5/897, 87/124,
                             176/877.5)
  nr_matching_genes_orig =  c( 121/133,175/177,174/177, 137/155.5,161/167, 153/196.5,
                               172.5/179)
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(c("BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),2)
  type = c(rep("Corrected counts",7),rep("Uncorrected Counts",7))
  full_rat = data.frame(nr_matching_genes,methods,type)
  
  
  p3 = ggplot(full_rat, aes( y=nr_matching_genes, x=as.factor(methods),fill=type)) + 
    geom_bar(alpha=0.9,stat='identity',position = "dodge") +
    geom_text(aes(label = sprintf("%.3f", nr_matching_genes), y= nr_matching_genes),  vjust = 1.5,
              position =position_dodge(width = 0.9))+
    
    theme(plot.title = element_text(size=20, face="bold", 
                                    margin = margin(10, 0, 10, 0)))+#ggtitle('Ratio of cells that change cluster')+
    xlab("Method") + ylab("Percentage of matching genes") + labs(fill = "Method") + scale_fill_brewer(palette="Dark2")
  
  p3 <- p3 + theme(
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
  diffexp_dat_neuro <- create_diffexp_dat_neuro()
  diffexp_dat_orig_neuro <- create_diffexp_dat_orig_neuro()
  
  p0_neuro = ggplot(diffexp_dat_neuro, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.title = element_blank()
    
  )
  p2_neuro = ggplot(diffexp_dat_orig_neuro, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
    axis.text.x = element_text(angle=45, vjust=.9, hjust=1),
    plot.margin = margin(48, 5, 14, 0, "pt"),
    legend.title = element_blank()
    
  )
  
  c1 = c(386,816,817.5,363.5,755,283,805)
  c2 = diffexp_dat_neuro[2:8,1]
  nr_matching_genes_full = c1/c2
  c1 = c(386,825,817.5,363.5,800,794,820)
  c2 = diffexp_dat_orig_neuro[2:8,1]
  nr_matching_genes_orig = c1/c2
  
  
  
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(c("BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),2)
  type = c(rep("Corrected counts",7),rep("Uncorrected Counts",7))
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

combined_de = ((p1 + p2)/p3) | ((p1_neuro + p2_neuro)/p3_neuro)
combined_de = combined_de + plot_layout(guides="collect") + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9),legend.position = "bottom")
combined_de
ggsave("plots/de_combined.png", width = 15, height = 9, units = "in")


run_plots_simul_pbmc<- function(){
  diffexp_simul_pbmc <- create_diffexp_simul_pbmc()
  diffexp_simul_pbmc_orig <- create_diffexp_simul_pbmc_orig()
  p0_simul_pbmc = ggplot(diffexp_simul_pbmc, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
  
  p2_simul_pbmc = ggplot(diffexp_simul_pbmc_orig, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
  c1 = c(114,145,156,126.5,158,85,161)
  c2 = diffexp_simul_pbmc[2:8,1]
  nr_matching_genes_full = c1/c2
  nr_matching_genes_full
  c1 = c(114,158,156,126.5,150,145,156)
  c2 = diffexp_simul_pbmc_orig[2:8,1]
  nr_matching_genes_orig = c1/c2
  nr_matching_genes_orig
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(c("BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),2)
  type = c(rep("Corrected counts",7),rep("Uncorrected Counts",7))
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
  diffexp_simul_neuro <- create_diffexp_simul_neuro()
  diffexp_simul_neuro_orig <- create_diffexp_simul_neuro_orig()
  p1_simul_neuro = ggplot(diffexp_simul_neuro, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
  
  #p1_simul_neuro<- p0_simul_neuro + scale_y_break(c(1100,1200),ticklabels =c(1200,1500),scales=c(0.2),space = 0.1)# +
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
  p2_simul_neuro = ggplot(diffexp_simul_neuro_orig, aes( y=nr_genes, x=methods,fill=data_names)) + 
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
  c1 = c(389,806,816.5,356,754.5,182,807)
  c2 = diffexp_simul_neuro[2:8,1]
  nr_matching_genes_full = c1/c2
  nr_matching_genes_full
  c1 = c(389,814,816,356,805,602,817)
  c2 = diffexp_simul_neuro_orig[2:8,1]
  nr_matching_genes_orig = c1/c2
  nr_matching_genes_orig
  nr_matching_genes = c(nr_matching_genes_full,nr_matching_genes_orig)
  methods = rep(c("BBKNN", "Combat","Harmony","LIGER","MNN","SCVI","Seurat"),2)
  type = c(rep("Corrected counts",7),rep("Uncorrected Counts",7))
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

combined_de = ((p1_simul_pbmc + p2_simul_pbmc)/p3_simul_pbmc) | ((p1_simul_neuro + p2_simul_neuro)/p3_simul_neuro)
combined_de = combined_de + plot_layout(guides="collect") + plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=9),legend.position = "bottom")
combined_de
ggsave("plots/de_simulated_combined.png", width = 15, height = 9, units = "in")
