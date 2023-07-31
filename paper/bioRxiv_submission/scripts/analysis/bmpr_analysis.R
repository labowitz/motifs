## Change these variables
filename = "bmpr" # how name is appended to files for pathway
pathway_index = 4 # Tgfb pathway index

## Directories for the files
source("./scripts/analysis/imports.R")

output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")

pathway_name =  param_list[[pathway_index]][[1]] # Tgfb pathway name
min_genes_pathway = param_list[[pathway_index]][[2]] # Tgfb min. number of genes expressed
min_expr_threshold = 0.2 # Tgfb minimum expression threshold for gene to be on

optimal_k_pathway = 20 # optimal pathway components, computed from z-score

## Make a plot of k_opt profiles for cell types with pathway "ON
pipe_run <- quickPipeline(master_seurat = master_seurat, 
                          k_final = optimal_k_pathway, 
                          silent_plot = F, 
                          which_pathway = pathway_name, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold, 
                          which_profiles = 'both', 
                          save_pdf = T, 
                          pdf_file = paste(output_dir, "_profiles.pdf", sep = ""))

pipe_run$df_devel %>% select(cell_id, class_label) %>% write.csv(paste("./scripts/analysis/outputs/", filename, "_profiles.csv", sep = ""))

## Computing the optimal silhouette score
silh_plt = silhouettePlot(which_pathway = pathway_name, 
                          min_ON = min_genes_pathway, 
                          min_expr = min_expr_threshold, 
                          n_bootstraps = 100)

g <- silhouette_zscore(silh_result = silh_plt,
                       pathway_name = pathway_name,
                       x_offset = 6,
                       max.y = 0.43,
                       min.y =0.17 # adjust axis parameters 
)
g
ggsave(paste(output_dir, "_silh_z-score.pdf", sep = ""), width = 7, height = 7)