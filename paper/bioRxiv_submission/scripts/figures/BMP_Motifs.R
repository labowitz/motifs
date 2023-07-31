## Change these variables
filename = "tgfb" # how name is appended to files for pathway
pathway_index = 2 # Tgfb pathway index

## Directories for the files
source("./scripts/analysis/imports.R")

output_dir = paste("./scripts/figures/", filename, sep = "")

pathway_name =  param_list[[pathway_index]][[1]] # Tgfb pathway name
min_genes_pathway = param_list[[pathway_index]][[2]] # Tgfb min. number of genes expressed
min_expr_threshold = 0.2 # Tgfb minimum expression threshold for gene to be on

optimal_k_pathway = 30 # optimal pathway components, computed from z-score
ptimal_k_pathway = 30 # optimal pathway components, computed from z-score

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
                 height = 10,
                 width = 5)
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
