library(zoo)
library(pracma)
library(dplyr)
library(ggplot2)
library(gridExtra)
## Directories for the files
library(pals)
library(tidyverse)

source("./scripts/analysis/imports.R")
# Load Pathbank analysis results
load("~/Documents/Elowitz Lab/pathway_motifs/data/processed_data/pathbank_silhouetteControl2022.rda")

filename = "eph-ephrin" # how name is appended to files for pathway
output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")

all_pathways$pathway[all_pathways$pathway=='Eph_r'] <- 'Eph-Ephrin'
all_pathways$pathway[all_pathways$pathway=='Eph_l'] <- 'Eph-Ephrin'

data_real_plot$pathway[data_real_plot$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
data_real_plot$pathway[data_real_plot$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
data_real_plot$pathway[data_real_plot$pathway=='Srsf'] <- 'RNA-splicing by SR protein family'
data_real_plot$pathway[data_real_plot$pathway=='Wnt'] <- 'Frizzled and Lrp5/6 receptors for Wnt/B-Catenin Signaling'
data_real_plot$pathway[data_real_plot$pathway=='Fgfr'] <- 'FGF cell signaling proteins'

silh_plt <- readRDS(paste(output_dir, "silh_plt.RDS", sep=""))

# Parameters to show
pathway_name = "Eph-Ephrin"

min_expr_threshold <- 0.3
min_genes_pathway <- 2

# Plot densities
plotDensity <- function(){
  # Make master DF and plot 
  all_data_pathbank <- rbind(data_real_rank %>% dplyr::select(pathway, n_genes, pathway_class, diff, max_z), 
                             data_random_diff %>% dplyr::select(pathway, n_genes, pathway_class, diff, max_z)) 
  
  all_data_pathbank[all_data_pathbank$pathway_class!='random',]$pathway_class <- 'pathway'
  
  
  p5 <- all_data_pathbank %>% ggplot(aes(x = diff, fill = pathway_class)) +geom_density(alpha = 0.6) + theme_minimal()  + theme(legend.position = "none") + 
    geom_vline(xintercept = quantile(data_random_diff$diff, 0.01)) + 
    coord_cartesian(xlim = c(-4,5))
  
  p6 <- data_real_rank %>% ggplot(aes(x = diff)) + geom_histogram() + theme_minimal() +
    geom_vline(xintercept = quantile(data_random_diff$diff, 0.01)) + 
    geom_vline(xintercept = quantile(data_random_diff$diff, 0.99)) + 
    coord_cartesian(xlim = c(-4,5))
  return(list(p5,p6))	
}

# Individual histograms
clusterScoreDist <- function(which_pathway = 'Bmp_Tgfb'){
  pathway_size <- all_pathways %>% dplyr::filter(pathway== which_pathway) %>% pull(gene) %>% length 
  null_dist <-  data_random_plot %>% dplyr::filter(n_genes == pathway_size)
  
  p1 <- null_dist %>% ggplot(aes(x = cluster_ratio)) + 
    geom_density(col = "blue") + 
    theme_minimal()
  
  path_stats <- data_real_plot %>% dplyr::filter(pathway==which_pathway) 
  
  p1 <- p1 + geom_vline(xintercept = path_stats$max_peak/path_stats$n_genes, col = 'blue') + xlab('N clusters / n genes')+ 
    theme(text = element_text(size = 10)) + 
    ggtitle(paste(which_pathway , ' vs ', dim(null_dist)[1],' random sets' , sep =""))
  return(p1) 
}

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
# 2. Take the mean 

mean_k_peak = median(z_peaks$max_k)
data_real_plot$pathway[data_real_plot$pathway=='Eph_l'] <- pathway_name
data_real_plot[data_real_plot$pathway==pathway_name,]$max_peak = maxk
data_real_plot[data_real_plot$pathway==pathway_name,]$n_genes = dim(all_pathways[all_pathways$pathway == pathway_name,])[[1]]
data_real_plot[data_real_plot$pathway==pathway_name,]$ratio = maxk/data_real_plot[data_real_plot$pathway==pathway_name,]$n_genes

data_real_plot <- data_real_plot %>% filter(pathway != "Eph_r")

p1 <- clusterScoreDist(pathway_name) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + theme(axis.line = element_line(size = 0.5),
                                                                                                                        panel.grid.minor.x = element_blank(),
                                                                                                                        panel.grid.minor.y = element_blank())

ggsave("./scripts/analysis/outputs/interpathway/null_dist_Eph-Ephrin.pdf")