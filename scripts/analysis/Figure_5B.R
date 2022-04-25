## Load in libraries
source("./scripts/analysis/imports.R")

library(zoo)
library(pracma)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(pals)
library(tidyverse)

## Load Pathbank analysis results
load("~/Documents/Elowitz Lab/pathway_motifs/data/processed_data/pathbank_silhouetteControl2022.rda")

## Directories for the outputs
output_dir = paste("./scripts/figures/", sep = "")

## Change the pathway names for display
all_pathways$pathway[all_pathways$pathway=='Eph_r'] <- 'Eph-Ephrin'
all_pathways$pathway[all_pathways$pathway=='Eph_l'] <- 'Eph-Ephrin'

data_real_plot$pathway[data_real_plot$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
data_real_plot$pathway[data_real_plot$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
data_real_plot$pathway[data_real_plot$pathway=='Srsf'] <- 'RNA-splicing by SR protein family'
data_real_plot$pathway[data_real_plot$pathway=='Wnt'] <- 'Frizzled and Lrp5/6 receptors for Wnt/B-Catenin Signaling'
data_real_plot$pathway[data_real_plot$pathway=='Fgfr'] <- 'FGF cell signaling proteins'
data_real_plot <- data_real_plot %>% filter(pathway != "Eph_r")

## Compute the whole combined Eph-ephrin pathway statistics (which were not originally computed)

# Set pathway parameters
pathway_name = "Eph-Ephrin"
min_expr_threshold <- 0.3
min_genes_pathway <- 2

# Read in the silhouette analysis results (to find the peak num of clusters)
silh_plt <- readRDS(paste0("./scripts/analysis/outputs/aws_pathways/", pathway_name, "_silh_plt.RDS", sep=""))

# Perform the computations for the Eph-ephrin pathway of 20 genes
n_peaks = 10
pct_peak = 0.7
k_smooth = 3

zz = silh_plt[[2]] %>% pull(m)

zz[is.na(zz)] = 0 

zz_smooth <- rollmean(zz, k = k_smooth)
zz[k_smooth:length(zz)] = zz_smooth # shift the array to compensate for offset 

z_peaks <- findpeaks(zz) %>% as.data.frame
names(z_peaks) <- c('peak','max_k','start','end')

z_peaks %>% arrange(desc(peak))  -> z_peaks 
z_peaks %>% top_n(n = n_peaks, wt = peak) -> z_peaks # select the top peaks 
# filter for peaks that deviate more than 20% from the maximum 
z_peaks %>% dplyr::filter(peak> pct_peak *max(peak)) -> z_peaks

# select the k with the max z-score from the peaks 
# 1. Take the k to the max z 
maxk = z_peaks %>% top_n(n = 1, wt = max_k) %>% pull(max_k)

data_real_plot$pathway[data_real_plot$pathway=='Eph_l'] <- pathway_name
data_real_plot[data_real_plot$pathway==pathway_name,]$max_peak = maxk
data_real_plot[data_real_plot$pathway==pathway_name,]$n_genes = dim(all_pathways[all_pathways$pathway == pathway_name,])[[1]]
data_real_plot[data_real_plot$pathway==pathway_name,]$ratio = maxk/data_real_plot[data_real_plot$pathway==pathway_name,]$n_genes

## Function that computes the deviation of the pathway from null model
deviationFromNull <- function(data_real_plot = data.frame() ,
                              summary_stats = data.frame() , 
                              use_smooth = F){
  
  data_real_plot %>% mutate(ratio = max_peak/n_genes) -> data_real_plot
  
  # Find smooth mean and standard deviation for each timepoint 
  
  diff_baseline = c() 
  k_window = 3
  summary_stats %>% 
    mutate(smooth_mean = rollmean(cluster_ratio,k=k_window,fill=NA), 
           smooth_lower_sd = rollmean(cluster_ratio-sd_k, k = k_window,fill=NA), 
           smooth_upper_sd = rollmean(cluster_ratio + sd_k, k=k_window,fill=NA)) -> smooth_stats
  
  # we can use smoothed mean and variance to compute the differences 
  # useful if the distributions are not large enough 	
  if(use_smooth){
    smooth_stats$cluster_ratio <- smooth_stats$smooth_mean # replace with smoothed version
    smooth_stats$sd_k          <- smooth_stats$smooth_mean - smooth_stats$smooth_lower_sd# replace with smoothed version
  }	
  
  # for each pathway in the data.frame find the corresponding null distribution 
  # and compute the difference to the mean (only 
  for(i in 1:dim(data_real_plot)[1]){
    pathway_size = data_real_plot[i,]$n_genes 
    if(pathway_size>7){
      baseline = smooth_stats[smooth_stats$n_genes==pathway_size,]$cluster_ratio
      baseline_sd = smooth_stats[smooth_stats$n_genes==pathway_size,]$sd_k
      
      diff_baseline[i] = (data_real_plot[i,]$ratio - baseline)/ baseline_sd
    }else{
      diff_baseline[i]= 0 
    }
  }
  # assign deviation from mean background 
  data_real_plot$diff <- diff_baseline
  
  return(data_real_plot) # should work for random sets and real pathways
}

## Function that computes the p-value of deviation compared to the null model
getPvalues <- function(data_real_plot  = data.frame() , 
                       data_random_plot  = data.frame() # control distribution
){
  # p-val based on current null distribution 
  p_vals = c() 
  for(i in 1:dim(data_real_plot)[1]){
    pathway_size <- data_real_plot[i,]$n_genes 
    pathway_clust_ratio <- data_real_plot[i,]$max_peak / data_real_plot[i,]$n_genes
    # get ratio distribution 
    aa <- data_random_plot %>% dplyr::filter(n_genes == pathway_size ) %>% pull(cluster_ratio)
    if(data_real_plot[i,]$diff <0 ){
      p_vals[i] <- sum(aa<pathway_clust_ratio)/length(aa)
    }else{
      p_vals[i] <- sum(aa>pathway_clust_ratio)/length(aa)
    }
  }
  
  p_vals[is.na(p_vals)] <- 1
  data_real_plot$pval <- p_vals
  ### DONE pval 
  
  return(data_real_plot)
}

## Function to make plot of pathway deviations from null model ranked
makeRankPlot <- function(data_real_plot= data.frame() ){
  # fix pathway names 
  
  data_real_plot$full_name <- paste(data_real_plot$pathway, " (", data_real_plot$n_genes,")", sep="")
  
  data_real_plot <- data_real_plot %>% 
    dplyr::filter(!pathway %in% c('Bmp','Bmp_fun',
                                  'Eph_l','Bmp_down') ) 
  data_real_rank <- data_real_plot %>% dplyr::filter(n_genes>6)
  
  p1 <- data_real_rank %>% mutate(full_name = fct_reorder(full_name , diff)) %>% 
    ggplot(aes(x = full_name , y = diff)) + 
    geom_point() + coord_flip() + theme_minimal()+ 
    theme(text = element_text(size = 10 ))
  
  data_sig <- data_real_rank %>% dplyr::filter(pval <=0.05)
  p1 <- p1 + geom_point(data = data_sig, aes(x = full_name , y = diff), color = 'blue') 
  return(p1) 
}

data_random_plot %>% group_by(n_genes) %>% summarise(cluster_ratio= mean(max_peak/n_genes), sd_k = sd(max_peak/n_genes)) -> summary_stats

deviationFromNull(data_real_plot , summary_stats = summary_stats, use_smooth = F) %>% 
  getPvalues(data_random_plot = data_random_plot) %>% makeRankPlot() -> dev_plot

ggsave(paste0(output_dir, "Figure_5B.pdf", sep=""))