## Change these variables
filename = "notch" # how name is appended to files for pathway
pathway_index = 1 # Notch pathway index

## Directories for the files
source("./scripts/analysis/imports.R")
output_dir = paste("./scripts/figures/", filename, sep = "")

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

# Diversity rank plots
rank_plot <- rank_diversity(
  which_pathway = pathway_name, 
  k_motifs = optimal_k_pathway, 
  min_expression = min_expr_threshold, 
  min_genes_pathway = min_genes_pathway,
  embedding_matrix = pca_proj
)
rank_plot
ggsave(paste(output_dir, "_rank_plot.pdf", sep = ""), width = 7, height = 7)

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

## Plot the most diverse profiles, aka the motifs
# 1. New: Select all profiles 
diverse_df <- control_res$profiles
# 2. tidy data.frame to average gene expression within a pathway profile 
diverse_df %>% pivot_longer(cols = genesPathway(pathway_name), 
                            names_to = 'gene', 
                            values_to = 'expression') %>% 
  select(cell_id, cell_ontology_class, Tissue, Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy 

# 3. average profile 
diverse_tidy %>% group_by(class_label,gene) %>% summarise(mean_expr = mean(expression), 
                                                          rank = mean(rank),diversity = mean(diversity), 
                                                          cell_types = mean(n)) -> diverse_tidy

# 4. wide format 
diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types), names_from=gene,values_from=mean_expr) %>% tibble::column_to_rownames('class_label')-> diverse_mat

# 5. We do the filtering here either for motifs or for non-diverse profiles 
control_res$diversity %>% 
  dplyr::filter(type=='transcriptome') %>% # choose the null model 
  pull(d) %>% quantile(0.90) -> divers_quantile

diverse_mat %>% dplyr::filter(diversity>divers_quantile) -> diverse_mat 

# annotation for significant pairs 
text_ann = matrix('',dim(diverse_mat[,genesPathway(pathway_name)])[1],dim(diverse_mat[,genesPathway(pathway_name)])[2])
text_ann[diverse_mat[,genesPathway(pathway_name)]>=0.1] = '+'

pdf(paste(output_dir, "_motif_profiles.pdf", sep = ""))
motif_heatmap <- superheat(diverse_mat[,genesPathway(pathway_name)],
                           pretty.order.rows = T,
                           X.text = text_ann,
                           X.text.size = 4,
                           X.text.col = "white",
                           #heat.pal = black_pal(100),
                           bottom.label.text.angle = 90, 
                           yr = sqrt(diverse_mat$cell_types),
                           yr.plot.type='bar',
                           yr.axis.name = "N cell types",
                           row.title = "Pathway motifs",
                           column.title = "Pathway genes",
                           bottom.label.size = 0.3,
                           grid.hline.col = "white",
                           grid.hline.size = 2, 
                           
)
dev.off()

# Get the developmental cell types in which Notch is expressed
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

devel_adult <- pipe_run$df_devel

devel_adult$cell_ontology_class %>% unique() %>% write.csv(paste(output_dir, "_all_cell_types.csv", sep = ""))

if(!'cell_id' %in% names(master_seurat@meta.data))
  devel_adult %>% mutate(cell_id = paste(global_cluster, dataset, sep = "_")) -> devel_adult

devel_adult$cell_ontology_class %>% unique() %>% write.csv(paste(output_dir, "_devel_cell_types.csv", sep = ""))

## quickPipeline processes the counts from an object with a specified optimal number of clusters
quickPipeline <- function(master_seurat = c(),              # Seurat object
                          which_pathway = 'Notch',          # Pathway name to draw from all_pathways object
                          k_final = 25,                     # k value for clusters
                          min_genes_on = 1,                 # Min. number of genes for pathway to be "ON"
                          min_expr = 0.2,                   # Min. expression cutoff for gene to be "ON"
                          which_profiles = 'devel',         # Include developmental, adult, or both datasets
                          rand_control = F,                 #
                          silent_plot = T,                  # Whether or not to return heatmap of data and class labels
                          manual_filter_cells = c(),        #
                          verbose = F,
                          save_pdf = F,
                          pdf_file = ""
){
  if(!rand_control){
    this_pathway = all_pathways %>% dplyr::filter(pathway == which_pathway) %>% pull(gene)
  }else{
    this_pathway = which_pathway # random set directly as input
  }
  
  # Get pathway genes in the dataset
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]
  
  # Get the MinMax scaled and normalized data and annotations
  df_devel <- normalizedDevel(master_seurat = master_seurat,
                              this_pathway,
                              sat_val = 0.99,
                              fill_zero_rows = F,
                              which_datasets = which_profiles)
  
  # Add a column of cell_type based on the different column values
  df_devel %>% mutate(cell_type = paste(global_cluster, '--',Tissue,': ', cell_ontology_class,'-', age, sep="")) -> df_devel
  
  # Compute the number of genes on in the pathway
  df_devel$genes_on = rowSums(df_devel[, this_pathway] > min_expr)
  
  row.names(df_devel) <- df_devel$cell_type
  
  # Two types of filtering.
  # Either by minimum expression of pathway genes OR user-specified list of cell types.
  if(length(manual_filter_cells) == 0){
    # Filter out cells with pathway "ON"
    df_devel %>% dplyr::filter(genes_on > min_genes_on) -> df_devel
  }else{
    # Else filter user-specified cell types
    df_devel %>% dplyr::filter(cell_id %in% manual_filter_cells) -> df_devel
    
    # minority of gene sets will have non-expressing cells
    # this happens very rarely and only 1 or 2 cells will have no expression, so we can remove them
    expressing_cells = df_devel$cell_id[rowSums(df_devel[, this_pathway]) > 0]
    df_devel %>% dplyr::filter(cell_id %in% expressing_cells) -> df_devel
  }
  
  # Heatmap to for computing the distance tree for cell types
  p = pheatmap(df_devel[,this_pathway ],
               #annotation_row = df_devel %>% select(dataset, Cell_class),
               #annotation_colors = colors_1206,
               show_rownames = T,
               fontsize = 5,
               cutree_rows = k_final,
               clustering_method = 'ward.D2',
               clustering_distance_rows = dist.cosine(as.matrix(df_devel[, this_pathway])),
               cluster_cols = F,
               silent = T)
  
  # Get class labels
  cos_labels = cutree(p$tree_row, k = k_final) # Cosine tree at specified k-cut
  df_devel$class_label = cos_labels # Map class labels
  
  df_devel$class_label <- as.character(df_devel$class_label)
  
  # Map colors to the class_labels and store in list called annoCols
  annoCols <- makeQualitativePal(k_final, glasbey_use = T, skip = 1) # Make Glasbey palette with k_final categorical colors
  names(annoCols) <- unique(df_devel$class_label)
  annoCols <- list("class_label" = annoCols)
  annoCols$Cell_class <- colors_1206$Cell_class
  annoCols$age <- colors_1206$age
  annoCols$dataset <- colors_1206$dataset
  
  # Heatmap to return with class_label annotations
  if(!save_pdf){
    p = pheatmap(df_devel[, this_pathway],
                 show_rownames = F,
                 fontsize = 5,
                 annotation_row = df_devel %>% select(class_label),
                 annotation_colors = annoCols,col = black_pal(100),
                 cutree_rows = k_final, clustering_method = 'ward.D2',
                 clustering_distance_rows = dist.cosine(as.matrix(df_devel[, this_pathway])),
                 cluster_cols = F,
                 silent = silent_plot)
  }else{
    p = pheatmap(df_devel[, this_pathway],
                 show_rownames = F,
                 fontsize = 5,
                 annotation_row = df_devel %>% select(class_label, Cell_class, age, dataset),
                 annotation_colors = annoCols,
                 col = black_pal(100),
                 cutree_rows = k_final,
                 clustering_method = 'ward.D2',
                 clustering_distance_rows = dist.cosine(as.matrix(df_devel[, this_pathway])),
                 cluster_cols = F,
                 silent = silent_plot,
                 filename = pdf_file,
                 height = 11,
                 width = 8.5)
  }
  
  df_devel$class_label <- as.numeric(df_devel$class_label)
  
  df_devel %>% gather(key = 'gene', value ='expression', this_pathway) %>%
    group_by(class_label, gene) %>%
    summarise(mean_expr = mean(expression), n_cell_types = n()) %>%
    spread(gene, mean_expr) %>% tibble::column_to_rownames(var = "class_label") -> x
  
  # Return count matrix, annotations, and heatmap for selected pathway
  return(list('matrix' = x, 'df_devel' = df_devel, 'profiles' = p))
}

pipe_run <- quickPipeline(master_seurat = master_seurat, 
                          k_final = optimal_k_pathway, 
                          silent_plot = F, 
                          which_pathway = pathway_name, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold, 
                          which_profiles = 'both', 
                          save_pdf = T, 
                          pdf_file = paste(output_dir, "_profiles.pdf", sep = ""))
