source("./scripts/analysis/imports_new.R")
library(networkD3)
library(webshot)
library(stringr)

## Directories for the files
fig_dir = "./scripts/figures/"
pathway_df = readRDS("./data/processed_data/all_pathways.RDS")

pathway_name =  "Tgf-beta family receptors" # tgfb pathway name
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj = master_seurat)

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
optimal_k_pathway = 30 # optimal pathway components, computed from z-score
diverse_quantile = 0.9

## Running fullControlPathway
control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                  k_final = optimal_k_pathway,
                                  seurat_obj = master_seurat, # seurat object
                                  null_list = hvg_genes, #list of highly variable genes 
                                  n_samples = 100, 
                                  filter_manual = T,
                                  min_genes_on = min_genes_pathway, 
                                  min_expr = min_expr_threshold, 
                                  n_pcs = 100, # how many PCs to use
                                  manual_embedding = pca_proj, # PCA embedding for diversity 
                                  dist_metric = "euclidean"
)

## Get profile labels of all cell types
glob_dendr <- global_dendr(control_res = control_res,
                           seurat_obj = master_seurat,
                           use_pca = T,
                           n_pcs = 1:20,
                           clust_method = 'ward.D2',
                           dist_metric ='cosine')

## Construct the dataframe for Sankey plot

# Generate colors for 31 pathway profiles (includes non-expressing profile 0)
colors_1206$class_label = makeQualitativePal(length(glob_dendr[[2]]$class_label %>% unique) + 1)
label_order <- c(21, 10, 24, 9, 8, 14, 13, 16, 7, 20, 15, 23, 22, 27)
glob_dendr[[2]] <- glob_dendr[[2]] %>% dplyr::filter(class_label %in% label_order)

glob_dendr[[2]]$class_label <- as.numeric(glob_dendr[[2]]$class_label)
temp <-rep(NA, length(glob_dendr[[2]]$class_label))
# Match label order of dendrogram in Fig. 3B
str_col = ""

## Printing string of colors correctly
for (i in 1:length(label_order)){
  temp[glob_dendr[[2]]$class_label == label_order[[i]]] <- i
  str_col <- cat(str_col, "\"", colors_1206$class_label[[sort(label_order)[[i]] + 1]], "\", ", sep="")
}

glob_dendr[[2]]$class_label <- temp

# Format dataframe for Sankey plot
# For each profile, count how many cell types are in each tissue class
cell_types = rep(c("Epithelium", "Macrophage", "Fibroblast", "Endothelial"), 
                 each = length(label_order))
receptor_profiles = rep(1:length(label_order), times = 4)
vals <- rep(0, length(label_order)*4)

# DataFrame that will store the counts
# Column cell types is the tissue profile
# Column receptor profile is the numerical receptor profile (0-30)
# Sankey takes only numerical values for nodes, so we rename the nodes
# 0-3 are the tissues, 4-34 are the profiles.
df <- data.frame(cell_types = cell_types, 
                 receptor_profiles = receptor_profiles, 
                 n = vals,
                 IDsource=rep(NA,length(label_order)*4),
                 IDtarget=rep(NA,length(label_order)*4))

# Count how many cell types are in Epithelium for each data frame
df_epi <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class == "Epithelium") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Epithelium" & df$receptor_profiles %in% df_epi$class_label,]$n <- df_epi$n
df[df$cell_types == "Epithelium",]$IDtarget <- 1:length(label_order) + 3
df[df$cell_types == "Epithelium",]$IDsource <- rep(0, length(label_order))

df_macro <- glob_dendr[[2]] %>% 
  dplyr::filter(str_detect(cell_ontology_class, fixed("macrophage", ignore_case = TRUE))) %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Macrophage" & df$receptor_profiles %in% df_macro$class_label,]$n <- df_macro$n
df[df$cell_types == "Macrophage",]$IDtarget <- 1:length(label_order) + 3
df[df$cell_types == "Macrophage",]$IDsource <- rep(1, length(label_order))

df_fibro <- glob_dendr[[2]] %>% 
  dplyr::filter(str_detect(cell_ontology_class, 
                           fixed("fibroblast", ignore_case = TRUE))) %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Fibroblast" & df$receptor_profiles %in% df_fibro$class_label,]$n <- df_fibro$n
df[df$cell_types == "Fibroblast",]$IDtarget <- 1:length(label_order) + 3
df[df$cell_types == "Fibroblast",]$IDsource <- rep(2, length(label_order))

df_endo <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class == "Endothelial") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Endothelial" & df$receptor_profiles %in% df_endo$class_label,]$n <- df_endo$n
df[df$cell_types == "Endothelial",]$IDtarget <- 1:length(label_order) + 3
df[df$cell_types == "Endothelial",]$IDsource <- rep(3, length(label_order))

names(df) = c("source", "target", "value", "IDsource", "IDtarget")

# Generate dataframe of node names
nodes <- data.frame("name" = c(c("Epithelium", 
                                 "Macrophage", 
                                 "Fibroblast", 
                                 "Endothelial"), 
                               as.character(label_order)))

# Sankey wants removal of zero counts
df <- df %>% dplyr::filter(value != 0)
df$IDsource <- as.numeric(df$IDsource)
df$IDtarget <- as.numeric(df$IDtarget)

# Specify matching between node labels and colors
color_spec <- 'd3.scaleOrdinal() .domain(["Epithelium", "Macrophage", "Fibroblast", "Endothelial", "7", "8", "9", "10", "13", "14", "15", "16", "20", "21", "22", "23", "24", "27"]) .range(["gray", "gray", "gray", "gray", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F", "#14F9FF"])'

s <- sankeyNetwork(Links = df, 
                   Nodes = nodes,
                   Source = "IDsource", 
                   Target = "IDtarget",
                   Value = "value", 
                   NodeID = "name",
                   fontSize= 12, 
                   nodeWidth = 30, 
                   iterations = 0,
                   colourScale = color_spec)

saveNetwork(s, paste0(fig_dir, "Figure_4G_motif.html", sep=""))
webshot::webshot(paste0(fig_dir, "Figure_4G_motif.html", sep=""), paste0(fig_dir, "Figure_4G_motif.pdf", sep=""), vwidth = 1000, vheight = 900)