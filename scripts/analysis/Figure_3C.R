library(networkD3)
library(webshot)
library(stringr)
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

cell_types = rep(c("Epithelium", "Macrophage", "Fibroblast", "Endothelial"), each = 31)
receptor_profiles = rep(1:31, times = 4)
vals <- rep(0, 124)

df <- data.frame(cell_types = cell_types, receptor_profiles = receptor_profiles, n = vals)

df_epi <- glob_dendr[[2]] %>% filter(Cell_class == "Epithelium") %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Epithelium" & df$receptor_profiles %in% df_epi$class_label,]$n <- df_epi$n
df[df$cell_types == "Epithelium",]$receptor_profiles <- 1:31 + 3
df[df$cell_types == "Epithelium",]$cell_types <- rep(0, 31)

df_macro <- glob_dendr[[2]] %>% filter(str_detect(cell_ontology_class, fixed("macrophage", ignore_case = TRUE))) %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Macrophage" & df$receptor_profiles %in% df_macro$class_label,]$n <- df_macro$n
df[df$cell_types == "Macrophage",]$receptor_profiles <- 1:31 + 3
df[df$cell_types == "Macrophage",]$cell_types <- rep(1, 31)

df_fibro <- glob_dendr[[2]] %>% filter(str_detect(cell_ontology_class, fixed("fibroblast", ignore_case = TRUE))) %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Fibroblast" & df$receptor_profiles %in% df_fibro$class_label,]$n <- df_fibro$n
df[df$cell_types == "Fibroblast",]$receptor_profiles <- 1:31 + 3
df[df$cell_types == "Fibroblast",]$cell_types <- rep(2, 31)

df_endo <- glob_dendr[[2]] %>% filter(Cell_class == "Endothelial") %>% group_by(class_label) %>% select(class_label) %>% count() %>% data.frame()
df[df$cell_types == "Endothelial" & df$receptor_profiles %in% df_endo$class_label,]$n <- df_endo$n
df[df$cell_types == "Endothelial",]$receptor_profiles <- 1:31 + 3
df[df$cell_types == "Endothelial",]$cell_types <- rep(3, 31)

names(df) = c("source", "target", "value")
nodes <- data.frame("name" = c(c("Epithelium", "Macrophage", "Fibroblast", "Endothelial"), as.character(label_order)))
#nodes <- data.frame(name = unique(c(df$source,df$target)),stringsAsFactors=FALSE)

df <- df %>% filter(value != 0)
df$source <- as.numeric(df$source)
df$target <- as.numeric(df$target)

colors_1206$class_label <- makeQualitativePal(optimal_k_pathway + 1, glasbey_use = T, skip = 0) 
names(colors_1206$class_label) <- 0:(optimal_k_pathway)

temp <- rep(NA, optimal_k_pathway + 1)

for (i in 1:length(label_order)){
  temp[names(colors_1206$class_label) == label_order[[i]]] <- i
}

names(colors_1206$class_label) <- temp
sorting <- sort(as.numeric(names(colors_1206$class_label)))
colors_1206$class_label <- colors_1206$class_label[as.character(sorting)]

colors_1206$class_label <- unname(colors_1206$class_label)
colors_1206$class_label <- c(rep("gray", 4), colors_1206$class_label)

color_scale <- 'd3.scaleOrdinal.range(["gray", "gray", "gray", "gray", "#FFFFFF", "#FF00B6", "#FF0000", "#0000FF", "#A10300", "#1F9698", "#201A01", "#720055", "#9A4D42", "#009FFF", "#783FC1", "#02AD24", "#B1CC71", "#00FFBE", "#00479E", "#14F9FF", "#C8FF00", "#886C00", "#00FF00", "#000033", "#858567", "#93D4FF", "#FFB79F", "#766C95", "#FE8F42", "#FFD300", "#005300", "#DC5E93", "#F1085C", "#FFACFD", "#DD00FF"])'

s <- sankeyNetwork(Links = df, Nodes = nodes,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   fontSize= 12, nodeWidth = 30, iterations = 0)

saveNetwork(s, paste0(fig_dir, "Figure_3C.html", sep=""))
webshot::webshot(paste0(fig_dir, "Figure_3C.html", sep=""), paste0(fig_dir, "Figure_3C.pdf", sep=""), vwidth = 1000, vheight = 900)