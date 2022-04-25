library(ggplot2)
library(viridis)

## Change these variables
filename = "tgfb" # how name is appended to files for pathway
pathway_index = 2 # Tgfb pathway index

## Directories for the files
source("./scripts/analysis/imports.R")
output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")
fig_dir = "./scripts/figures/"

pathway_name =  param_list[[pathway_index]][[1]] # Tgfb pathway name
min_genes_pathway = param_list[[pathway_index]][[2]] # Tgfb min. number of genes expressed
min_expr_threshold = 0.2 # Tgfb minimum expression threshold for gene to be on

optimal_k_pathway = 30 # optimal pathway components, computed from z-score

## Running fullControlPathway
control_res  = fullControlPathway(this_pathway = pathway_name,
                                  k_pathway = optimal_k_pathway , 
                                  filter_pathway = 'both', # adult + devel 
                                  this_seurat = master_seurat, # seurat object
                                  null_list = hvg_genes, #list of highly variable genes 
                                  n_samples = 100, filter_manual = T,
                                  min_expr_gene = min_expr_threshold, 
                                  min_genes_ON = min_genes_pathway, # default > 2 genes ON 
                                  n_pcs_global = 100, # how many PCs to use
                                  embedding_matrix = pca_proj # PCA embedding for diversity 
)

# Plot global dendrogram with pathway states
glob_dendr <- globalDendrogramPlot2(
  control_res = control_res, 
  seurat_obj = master_seurat, 
  use_pca = T, 
  npcs = 1:20,
  clust_method = 'ward.D2',
  dist_metric ='cosine',
  save_dir = output_dir
)

glob_dendr[[2]]$class_label <- as.numeric(glob_dendr[[2]]$class_label)
temp <-rep(NA, length(glob_dendr[[2]]$class_label))
label_order <- c(0, 5, 2, 1, 26, 12, 18, 19, 9, 8, 11, 21, 14, 10, 28, 27, 22, 23, 3, 4, 25, 30, 24, 20, 16, 7, 6, 29, 15, 13, 17)

for (i in 1:length(label_order)){
  temp[glob_dendr[[2]]$class_label == label_order[[i]]] <- i
}

glob_dendr[[2]]$class_label <- temp

cell_types = rep(c("Epithelium", "Macrophage", "Fibroblast", "Endothelial", "Other"), each = 31)
receptor_profiles = rep(1:31, times = 5)
vals <- rep(0, 31*5)

df <- data.frame(cell_types = cell_types, receptor_profiles = receptor_profiles, n = vals)

df_epi <- glob_dendr[[2]] %>% filter(Cell_class == "Epithelium") %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Epithelium" & df$receptor_profiles %in% df_epi$class_label,]$n <- df_epi$n

df_macro <- glob_dendr[[2]] %>% filter(str_detect(cell_ontology_class, fixed("macrophage", ignore_case = TRUE))) %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Macrophage" & df$receptor_profiles %in% df_macro$class_label,]$n <- df_macro$n

df_fibro <- glob_dendr[[2]] %>% filter(str_detect(cell_ontology_class, fixed("fibroblast", ignore_case = TRUE))) %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Fibroblast" & df$receptor_profiles %in% df_fibro$class_label,]$n <- df_fibro$n

df_endo <- glob_dendr[[2]] %>% filter(Cell_class == "Endothelial") %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Endothelial" & df$receptor_profiles %in% df_endo$class_label,]$n <- df_endo$n

df_other <- glob_dendr[[2]] %>% filter(Cell_class != "Epithelium") %>% filter(!str_detect(cell_ontology_class, fixed("macrophage", ignore_case = TRUE))) %>% filter(!str_detect(cell_ontology_class, fixed("fibroblast", ignore_case = TRUE))) %>% filter(Cell_class != "Endothelial")  %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Other" & df$receptor_profiles %in% df_other$class_label,]$n <- df_other$n

# Subtract 1 from the class labels
df$receptor_profiles = df$receptor_profiles - 1

p <- ggplot(df, aes(fill=cell_types, y=n, x=receptor_profiles)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("Epithelium"=colors_1206$Cell_class["Epithelium"][[1]], "Macrophage"=colors_1206$Cell_class["Blood"][[1]], "Fibroblast"=colors_1206$Cell_class["Connective"][[1]], "Endothelial"=colors_1206$Cell_class["Endothelial"][[1]], "Other"="#D3D3D3")) +
  ggtitle("Profile Tissue Composition") +
  xlab("Profile Number") + ylab("Percent Composition") + 
  theme_bw()

ggsave(paste0(fig_dir, "Figure_3C_Inset.pdf", sep=""))