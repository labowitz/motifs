## Module with functions
source("./module/module_new.R")

## Seurat object to import
master_seurat <- readRDS("./data/processed_data/master_seurat.RDS") # Most recent Seurat object
pca_proj <- readRDS("./data/processed_data/pca_proj.RDS") # PCA projection from 1928 HVGs
hvg_genes <- readRDS("./data/processed_data/hvg_genes.RDS") # list of HVGs (shouldn't be used)

## Set of files to import 
pathway_list <- readRDS("./data/processed_data/pathway_list.RDS") # List of pathway components
param_list <- readRDS("./data/processed_data/param_list.RDS") # List of pathway parameters
all_pathways <- readRDS("./data/processed_data/all_pathways.RDS") # List of all pathways

## Colors to import
colors_1206 <- readRDS("./data/processed_data/colors_1206.RDS") # List with labels and corresponding colors
glasbey <- readRDS("./data/processed_data/glasbey.RDS") # List of glasbey colors for categorical labels