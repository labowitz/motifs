## Module with functions
source("./module/module.R")

## Seurat object to import
master_seurat <- readRDS("./data/processed_data/master_seurat.RDS") # Most recent Seurat object
pca_proj <- readRDS("./data/processed_data/umap_may.RDS") # PCA projection from 1928 HVGs
hvg_genes <- readRDS("./data/processed_data/hvg_genes.RDS") # list of HVGs (shouldn't be used)


## Set of files to import 
pathway_list <- readRDS("./data/processed_data/pathway_list.RDS") # List of pathway components
param_list <- readRDS("./data/processed_data/param_list.RDS") # List of pathway parameters
all_pathways <- readRDS("./data/processed_data/all_pathways.RDS") # List of all pathways


## Colors to import
colors_1206 <- readRDS("./data/processed_data/colors_1206.RDS") # List with labels and corresponding colors
glasbey <- readRDS("./data/processed_data/glasbey.RDS") # List of glasbey colors for categorical labels


# Alter the metadata of the file to reflect changes in "Organ_specific" labels
organ_specific <- read.table('./data/raw_data/organ_specific_metadata - organ_specific_metadata.tsv', header = T, sep = "\t")

# Save previous 
master_seurat@meta.data$Cell_class_prev <- master_seurat@meta.data$Cell_class

# Check indexes are the same 
#identical(row.names(master_seurat@meta.data[ master_seurat@meta.data$Cell_class=='Organ_specific',  ]) , organ_specific$index ) 

# replace with new Cell type classes 
master_seurat@meta.data[master_seurat@meta.data$Cell_class == 'Organ_specific', ]$Cell_class <- organ_specific$Cell_class