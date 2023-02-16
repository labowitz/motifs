source("./scripts/analysis/imports_new.R")
library(networkD3)
library(webshot)
library(stringr)

## Directories for the files
fig_dir = "./scripts/figures/"
pathway_df = readRDS("./data/processed_data/all_pathways.RDS")

pathway_name =  "Bmp_Tgfb" # tgfb pathway name
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df)

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
glob_dendr[[2]]$class_label <- as.numeric(glob_dendr[[2]]$class_label)
temp <-rep(NA, length(glob_dendr[[2]]$class_label))
# Match label order of dendrogram in Fig. 3B
label_order <- c(0, 5, 2, 1, 26, 12, 18, 19, 9, 8, 11, 21, 14, 10, 28, 27, 22, 23, 3, 4, 25, 30, 24, 20, 16, 7, 6, 29, 15, 13, 17)
str_col = ""

## Printing string of colors correctly
for (i in 1:length(label_order)){
  temp[glob_dendr[[2]]$class_label == label_order[[i]]] <- i
  str_col <- cat(str_col, "\"", colors_1206$class_label[[i]], "\", ", sep="")
}

glob_dendr[[2]]$class_label <- temp

# Format dataframe for Sankey plot
# For each profile, count how many cell types are in each tissue class
cell_types = rep(c("Epithelium", "Macrophage", "Fibroblast", "Endothelial"), 
                 each = 31)
receptor_profiles = rep(1:31, times = 4)
vals <- rep(0, 124)

# DataFrame that will store the counts
# Column cell types is the tissue profile
# Column receptor profile is the numerical receptor profile (0-30)
# Sankey takes only numerical values for nodes, so we rename the nodes
# 0-3 are the tissues, 4-34 are the profiles.
df <- data.frame(cell_types = cell_types, 
                 receptor_profiles = receptor_profiles, 
                 n = vals,
                 IDsource=rep(NA,31*4),
                 IDtarget=rep(NA,31*4))

# Count how many cell types are in Epithelium for each data frame
df_epi <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class == "Epithelium") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Epithelium" & df$receptor_profiles %in% df_epi$class_label,]$n <- df_epi$n
df[df$cell_types == "Epithelium",]$IDtarget <- 1:31 + 3
df[df$cell_types == "Epithelium",]$IDsource <- rep(0, 31)

df_macro <- glob_dendr[[2]] %>% 
  dplyr::filter(str_detect(cell_ontology_class, fixed("macrophage", ignore_case = TRUE))) %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Macrophage" & df$receptor_profiles %in% df_macro$class_label,]$n <- df_macro$n
df[df$cell_types == "Macrophage",]$IDtarget <- 1:31 + 3
df[df$cell_types == "Macrophage",]$IDsource <- rep(1, 31)

df_fibro <- glob_dendr[[2]] %>% 
  dplyr::filter(str_detect(cell_ontology_class, 
                    fixed("fibroblast", ignore_case = TRUE))) %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Fibroblast" & df$receptor_profiles %in% df_fibro$class_label,]$n <- df_fibro$n
df[df$cell_types == "Fibroblast",]$IDtarget <- 1:31 + 3
df[df$cell_types == "Fibroblast",]$IDsource <- rep(2, 31)

df_endo <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class == "Endothelial") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()
df[df$cell_types == "Endothelial" & df$receptor_profiles %in% df_endo$class_label,]$n <- df_endo$n
df[df$cell_types == "Endothelial",]$IDtarget <- 1:31 + 3
df[df$cell_types == "Endothelial",]$IDsource <- rep(3, 31)

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
color_spec <- 'd3.scaleOrdinal() .domain(["Epithelium", "Macrophage", "Fibroblast", "Endothelial", "0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30"]) .range(["gray", "gray", "gray", "gray", "#FFFFFF", "#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300", "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#DD00FF", "#201A01", "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F", "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF"])'

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

saveNetwork(s, paste0(fig_dir, "Figure_3C.html", sep=""))
webshot::webshot(paste0(fig_dir, "Figure_3C.html", sep=""), paste0(fig_dir, "Figure_3C.pdf", sep=""), vwidth = 1000, vheight = 900)