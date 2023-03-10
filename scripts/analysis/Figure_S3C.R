## Directories for the files
source("./scripts/analysis/imports_new.R")
fig_dir = "./scripts/figures/"
pathway_df = readRDS("./data/processed_data/all_pathways.RDS")

pathway_name =  "Tgf-beta family receptors" # tgfb pathway name
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj=master_seurat)

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

glob_dendr[[2]]$class_label <- as.numeric(glob_dendr[[2]]$class_label)

cell_types = rep(c("Epithelium", 
                   "Macrophage", 
                   "Fibroblast", 
                   "Endothelial", 
                   "Other"), 
                 each = 31)
receptor_profiles = rep(0:30, times = 5)
vals <- rep(0, 31*5)

df <- data.frame(cell_types = cell_types, 
                 receptor_profiles = receptor_profiles, 
                 n = vals)

df_other <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class != "Epithelium") %>% 
  dplyr::filter(!str_detect(cell_ontology_class, 
                            fixed("macrophage", 
                                  ignore_case = TRUE))) %>% 
  dplyr::filter(!str_detect(cell_ontology_class, 
                            fixed("fibroblast", 
                                  ignore_case = TRUE))) %>% 
  dplyr::filter(Cell_class != "Endothelial") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()

df[df$cell_types == "Other" & df$receptor_profiles %in% df_other$class_label,]$n <- df_other$n

df_epi <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class == "Epithelium") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()

df[df$cell_types == "Epithelium" & df$receptor_profiles %in% df_epi$class_label,]$n <- df_epi$n

df_macro <- glob_dendr[[2]] %>% 
  dplyr::filter(str_detect(cell_ontology_class, 
                    fixed("macrophage", 
                          ignore_case = TRUE))) %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()

df[df$cell_types == "Macrophage" & df$receptor_profiles %in% df_macro$class_label,]$n <- df_macro$n

df_fibro <- glob_dendr[[2]] %>% 
  dplyr::filter(str_detect(cell_ontology_class, 
                    fixed("fibroblast", 
                          ignore_case = TRUE))) %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()

df[df$cell_types == "Fibroblast" & df$receptor_profiles %in% df_fibro$class_label,]$n <- df_fibro$n

df_endo <- glob_dendr[[2]] %>% 
  dplyr::filter(Cell_class == "Endothelial") %>% 
  dplyr::group_by(class_label) %>% 
  dplyr::select(class_label) %>% 
  dplyr::count() %>% 
  data.frame()

df[df$cell_types == "Endothelial" & df$receptor_profiles %in% df_endo$class_label,]$n <- df_endo$n

# Subtract 1 from the class labels
df$receptor_profiles = df$receptor_profiles - 1

p <- ggplot(df, 
            aes(fill=cell_types, y=n, x=receptor_profiles)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("Epithelium"=colors_1206$Cell_class["Epithelium"][[1]], 
                             "Macrophage"=colors_1206$Cell_class["Blood"][[1]], 
                             "Fibroblast"=colors_1206$Cell_class["Connective"][[1]], 
                             "Endothelial"=colors_1206$Cell_class["Endothelial"][[1]], 
                             "Other"="#D3D3D3")) +
  ggtitle("Profile Tissue Composition") +
  xlab("Profile Number") + 
  ylab("Percent Composition") + 
  theme_bw()

ggsave(paste0(fig_dir, "Figure_3C_Inset.pdf", sep=""))