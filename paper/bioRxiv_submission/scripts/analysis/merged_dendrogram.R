## Directories for the files
filename = "tgfb_notch_wnt.pdf"
source("./scripts/analysis/imports.R")
output_dir = paste("./scripts/analysis/outputs/", filename, sep = "")

dist_metric = "cosine"

# Get the Tgfb pathway states
pathway_index = 2 # Tgfb pathway index
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

glob_dendr <- globalDendrogramPlot2(control_res = control_res, 
                                    seurat_obj = master_seurat, 
                                    use_pca = T, 
                                    npcs = 10, 
                                    save_pdf = F
)

df <- glob_dendr[[2]]
df <- df %>% select(-class_label)
df$tgfb_class <- glob_dendr[[2]]$class_label

# Get the Notch pathway states
pathway_index = 1 # Notch pathway index
pathway_name =  param_list[[pathway_index]][[1]] # Notch pathway name
min_genes_pathway = param_list[[pathway_index]][[2]] # Notch min. number of genes expressed
min_expr_threshold = 0.2 # Notch minimum expression threshold for gene to be on

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

glob_dendr <- globalDendrogramPlot2(control_res = control_res, 
                                    seurat_obj = master_seurat, 
                                    use_pca = T, 
                                    npcs = 10, 
                                    save_pdf = F
)

df$notch_class <- glob_dendr[[2]]$class_label

# Get the Wnt pathway states
pathway_index = 3 # Wnt pathway index
pathway_name =  param_list[[pathway_index]][[1]] # Wnt pathway name
min_genes_pathway = param_list[[pathway_index]][[2]] # Wnt min. number of genes expressed
min_expr_threshold = 0.2 # Wnt minimum expression threshold for gene to be on

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

glob_dendr <- globalDendrogramPlot2(control_res = control_res, 
                                    seurat_obj = master_seurat, 
                                    use_pca = T, 
                                    npcs = 10, 
                                    save_pdf = F
)

df$wnt_class <- glob_dendr[[2]]$class_label

colors_1206$class_label <- makeQualitativePal(optimal_k_pathway + 1, glasbey_use = T, skip = 0) # skip white color
names(colors_1206$class_label) <- 0:(optimal_k_pathway) %>% as.character()

colors_1206$tgfb_class <- colors_1206$class_label
colors_1206$wnt_class <- colors_1206$class_label
colors_1206$notch_class <- colors_1206$class_label

x_scaled <- glob_dendr[[1]]

if(dist_metric=='cosine') clust_dist_metric=dist.cosine(x_scaled) else clust_dist_metric = dist(x_scaled)

pheatmap(x_scaled, 
         annotation_row = df %>% select(tgfb_class, notch_class, wnt_class, Cell_class,age, dataset),
         annotation_colors = colors_1206, 
         show_colnames = F, 
         show_rownames = F,
         clustering_distance_rows = clust_dist_metric, 
         treeheight_col = 0,
         cutree_rows = 12 , 
         fontsize = 10,
         color = magma(100),
         filename = paste(output_dir,'_global_dendrogram.pdf',sep=""),
         height = 20, 
         width = 10
         )
