# Libraries for this analysis
library(infotheo)
library(aricode)
library(RColorBrewer)

# Directories for the files
source("./scripts/analysis/imports.R")
output_dir = "./scripts/analysis/outputs/interpathway/"

# Blues palette
blues_pal <- function(n){
  colfunc <- colorRampPalette(c("blue", "white"))
  return(colfunc(n))
}


# Import data
tgfb_niv <- read.csv("./scripts/analysis/outputs/tgfb_profiles.csv", header = T)
bmpr_niv <- read.csv("./scripts/analysis/outputs/bmpr_profiles.csv", header = T)
notch_niv <- read.csv("./scripts/analysis/outputs/notch_profiles.csv", header = T)
wnt_niv <- read.csv("./scripts/analysis/outputs/wnt_profiles.csv", header = T)
eph_ephr_niv <- read.csv("./scripts/analysis/outputs/eph-ephrin_profiles.csv", header = T)
eph_ab_niv <- read.csv("./scripts/analysis/outputs/eph_ab.csv", header = T)
rna_niv <- read.csv("./scripts/analysis/outputs/rna_splicing_sr.csv", header = T)
ubi_niv <- read.csv("./scripts/analysis/outputs/ubi_prot.csv", header = T)


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

all_cell_type$eph_ephr <- 0
all_cell_type[eph_ephr_niv$cell_id,]$eph_ephr <- eph_ephr_niv$class_label

all_cell_type$eph_ab <- 0
all_cell_type[eph_ab_niv$cell_id,]$eph_ab <- eph_ab_niv$class_label

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
               filename = paste("./scripts/analysis/outputs/interpathway/rel_freq_", name_i, "_", name_j, ".pdf", sep = ""))
      
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


pathwayConfusion(all_cell_type %>% select(tgfb, bmpr, notch, wnt, eph_ephr, ubi, rna), filename = "./scripts/analysis/outputs/interpathway/MI_ubi_rna_notch_eph_bmp_wnt.pdf", save = T)