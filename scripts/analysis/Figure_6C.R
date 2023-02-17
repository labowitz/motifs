# Libraries for this analysis
library(infotheo)
library(aricode)
library(RColorBrewer)

# Directories for the files
source("./scripts/analysis/imports_new.R")

data_dir = "scripts/figures/peak_analysis/silhouette_res/silh_rds/"
output_dir = "./scripts/figures/"

# Blues palette
blues_pal <- function(n){
  colfunc <- colorRampPalette(c("blue", "white"))
  return(colfunc(n))
}

# Import data
tgfb_niv <- read.csv(paste(data_dir, 'Tgf-beta family receptors_profiles.csv', sep=""))
bmpr_niv <- read.csv(paste(data_dir, "Bmp", "_profiles.csv", sep=""))
notch_niv <- read.csv(paste(data_dir, 
                   'Notch receptors, Dll ligands and Fringe proteins',
                   "_profiles.csv", sep=""))
wnt_niv <- read.csv(paste(data_dir, "Frizzled and Lrp5 6 receptors for Wnt B Catenin Signaling", 
                 "_profiles.csv", sep=""))
rna_niv <- read.csv(paste(data_dir, 'RNA-splicing by SR protein family', 
                 "_profiles.csv", sep=""))
ubi_niv <- read.csv(paste(data_dir, "Ubiquitin-Proteasome Pathway",
                 "_profiles.csv", sep=""))

# Add cell class labels to the master_seurat objet
all_cell_type <- master_seurat@meta.data %>% select(cell_id, dataset, age, cell_ontology_class)

all_cell_type$tgfb <- 0
all_cell_type[tgfb_niv$cell_id,]$tgfb <- tgfb_niv$class_label

all_cell_type$bmpr <- 0
all_cell_type[bmpr_niv$cell_id,]$bmpr <- bmpr_niv$class_label

all_cell_type$notch <- 0
all_cell_type[notch_niv$cell_id,]$notch <- notch_niv$class_label

all_cell_type$wnt <- 0
all_cell_type[wnt_niv$cell_id,]$wnt <- wnt_niv$class_label

all_cell_type$rna <- 0
all_cell_type[rna_niv$cell_id,]$rna <- rna_niv$class_label

all_cell_type$ubi <- 0
all_cell_type[ubi_niv$cell_id,]$ubi <- ubi_niv$class_label

# Interpathway correlations
pathwayConfusion <- function(motif_labels = data.frame(), 
                             clust_method = 'ward.D2',
                             save = F,
                             filename = ""){
  
  # We are testing different metrics for clustering comparison: 
  MI = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
  N_MI = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
  ari = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
  Adj_MI = matrix(0, dim(motif_labels)[2], dim(motif_labels)[2])
  
  for(i in 1:(dim(motif_labels)[2]-1)){
    for(j in (i+1):dim(motif_labels)[2]){
      
      name_i <- colnames(motif_labels)[[i]]
      name_j <- colnames(motif_labels)[[j]]
      
      long_pair_freq <- data.frame(all_cell_type %>% select(name_i, name_j) %>% group_by_at(c(name_i,name_j)) %>% count())
      
      mtx_pair_freq <- matrix(0, nrow = length(unique(long_pair_freq[,name_i])), ncol = length(unique(long_pair_freq[,name_j])))

      rownames(mtx_pair_freq) <- sort(unique(long_pair_freq[,name_i]))
      colnames(mtx_pair_freq) <- sort(unique(long_pair_freq[,name_j]))
      
      for(idx in 1:nrow(long_pair_freq)){
        row = long_pair_freq[idx, name_i]
        col = long_pair_freq[idx, name_j]
        
        mtx_pair_freq[toString(row), toString(col)] <- long_pair_freq[idx,"n"]
      }
      
      pheatmap(mtx_pair_freq, 
               cluster_rows = F, 
               cluster_cols = F,
               main = paste(name_i, " and ", name_j, sep = ""),
               filename = paste("data/processed_data/interpathway_corr/rel_freq_", name_i, "_", name_j, ".pdf", sep = ""))
      
      # Mutual information 
      MI[i,j] = infotheo::mutinformation(motif_labels[,i], motif_labels[,j])
      # Adjusted mutual information 
      Adj_MI[i,j] = AMI(motif_labels[,i], motif_labels[,j])
      # Normalized mutual information 
      N_MI[i,j] = NMI(motif_labels[,i], motif_labels[,j])
      # Adjusted random index 
      #ari[k] = adjustedRandIndex(motif_labels[,i], motif_labels[,j]) # Adjusted Rand Index
    }
  }
  resMI = Adj_MI
  resMI  = resMI + t(resMI)
  diag(resMI) <- 1
  resMI[lower.tri(resMI)] <- NA
  
  row.names(resMI)<-colnames(motif_labels)
  colnames(resMI)<-colnames(motif_labels)
  if(!save){
    pheatmap(resMI, 
             fontsize = 14, 
             col = magma(100), 
             cluster_cols = F, 
             cluster_rows = F, 
             na_col="white"
             )
  }
  else{
    pheatmap(resMI, 
             fontsize = 14, 
             col = magma(100), 
             filename = filename,
             cluster_cols = F, 
             cluster_rows = F
             )
    
  }
  return(resMI)    
}

# Plot the interpathway correlations
pathwayConfusion(all_cell_type %>% select(tgfb, bmpr, notch, wnt, ubi, rna), 
                 filename = paste0(output_dir, "Figure_6C.pdf", sep=""), 
                 save = T)