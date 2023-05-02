source("./scripts/analysis/imports.R")

# Specify pathway name (as saved in pathway_df) and get genes.
pathway_name =  'Tgf-beta family receptors'
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj = master_seurat)

# Thresholds for processing data.
min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
diverse_quantile = 0.9
optimal_k_pathway = 30

# Computes the dispersion on cell types expressing pathway using distance in PCA space
control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                  k_final = optimal_k_pathway,
                                  seurat_obj = master_seurat,
                                  n_samples = 100,
                                  filter_manual = T,
                                  min_genes_on = min_genes_pathway, 
                                  min_expr = min_expr_threshold, 
                                  n_pcs = 100,
                                  manual_embedding = pca_proj, # PCA embedding for dispersion 
                                  dist_metric = "euclidean"
)

div_res_mean_pairwise <- diverseFilt(control_res = control_res,
                       pathway_genes = pathway_genes,
                       type="motif"
)$diverse_mat %>% rownames()

## MAX pairwise distances
globalClustering <- function(df_devel = data.frame(),   # data.frame in wide format with genes and class labels
                             pathway_genes = c(),       # List of pathway genes
                             seurat_obj = c(),          # Seurat object
                             k_final = 25,              # k value for clusters
                             n_pcs = 100,               # Number of PCs to use
                             manual_embedding = c(),    # Matrix should represent a latent space (i.e. PCA, UMAP), do not provide if using n_pcs
                             dist_metric = 'euclidean'  # Distance metric
){
  
  df_devel$class_label <- df_devel$class_label %>% as.character()
  
  # Rename df_devel by cell_id
  row.names(df_devel) <- df_devel$cell_id
  
  # Set class_label colors
  colors_1206$class_label = makeQualitativePal(length(df_devel$class_label %>% unique))
  names(colors_1206$class_label) <- df_devel$class_label %>% unique
  
  # Get the Embedding
  # umap_coords by default uses cell_id as the row.name
  if(length(manual_embedding) == 0){
    # If user specifies n_pcs
    umap_coords <- Embeddings(seurat_obj, reduction = 'pca')
    umap_coords <- umap_coords[, 1:n_pcs]
    
  }else{
    # Else use the user's Embedding, when provided directly, n_pcs is ignored
    # row.names equal cell_id
    umap_coords = manual_embedding
    
  }
  
  # Use the cells that only are in our data.frame object, with the pathway "ON"
  scaled_data = umap_coords[df_devel$cell_id, ]
  
  # Compute distance in embedding space between cells
  if(dist_metric =='euclidean')  user_dist = dist(scaled_data)  else  user_dist = dist.cosine(scaled_data)
  
  # Global clustering of cell types uses euclidean distance + ward.D2
  p_global = pheatmap(scaled_data %>% t ,
                      annotation_col = df_devel %>% 
                        dplyr::select(age, dataset, Tissue, Cell_class, class_label),
                      annotation_colors = colors_1206,
                      clustering_method = 'ward.D2',
                      cutree_cols = k_final,
                      show_rownames = F,
                      show_colnames = F,
                      fontsize = 12,
                      silent = T,
                      clustering_distance_cols = user_dist)
  
  # For each pathway label, we calculate the mean distance in the embedding space for the cell types in that label
  
  k_final = length(df_devel$class_label %>% unique)
  
  # Distance metric to be used as matrix
  # We convert the distance to a matrix so we can compute the pair-wise similarity
  dist_mat = as.matrix(user_dist)
  
  # For each class, compute the pairwise similarity
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$class_label == x]
    
    if(length(class_cells) > 1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% max())
      
    }else{      # If there is only one member in the class, then pairwise similarity is zero
      return(0)
    }
    
  }) -> diversity_pathway
  
  # Perform the same pairwise similarity computation on global labels
  
  global_labels = cutree(p_global$tree_col, 
                         k = k_final)
  df_devel$global_label = global_labels
  
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$global_label == x]
    
    if(length(class_cells)>1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% max())
      
    }else{
      return(0)
    }
    
  }) -> diversity_global
  
  # Transcriptome diversity
  df_devel %>% 
    group_by(global_label) %>% 
    count %>% 
    as.data.frame() -> global_stats
  global_stats$diversity = diversity_global
  
  # Pathway diversity
  df_devel %>% 
    group_by(class_label) %>% 
    count %>% 
    as.data.frame() %>%
    mutate(class_label = as.numeric(class_label)) %>%
    arrange(class_label) -> path_stats
  path_stats$diversity = diversity_pathway
  
  # Compute silhouette scores for this clustering
  # By default, we perform gene clustering using the cosine distance and ward.D2
  # so we also use cosine distance for the silhouette score
  
  pathway_silh = clusteringSilhouette(df = df_devel, 
                                      dist_metric = 'cosine', 
                                      pathway_genes = pathway_genes)
  
  # Join with the main data.frame
  path_stats %>% left_join(pathway_silh, by = 'class_label') -> path_stats
  # Dummy column to match later
  global_stats$mean_silh = 0
  
  path_stats$type = "pathway"
  global_stats$type = "transcriptome"
  
  path_stats %>% rename(label = class_label) -> path_stats
  global_stats %>% rename(label = global_label) -> global_stats
  
  # Rank the profiles by diversity
  path_stats$rank = rank(path_stats$diversity)
  
  return(list('pathway' = path_stats, 'global' = global_stats))
}
fullControlPathway <- function(pathway_genes = c(),      # List of pathway genes
                               k_final = c(),            # k value for clusters
                               seurat_obj = c(),         # Seurat object
                               n_samples = 100,          # No. of samples to estimate positive control
                               filter_manual = F,
                               min_expr = 0.2,           # Min. expression cutoff for gene to be "ON"
                               min_genes_on = 1,         # Min. number of genes for pathway to be "ON"
                               n_pcs = 100,              # Number of PCs to use
                               manual_embedding = c(),   # Number of principal components to consider for diversity metric
                               dist_metric = 'euclidean' # distance metric to be use in the global embedding (e.g., PCA) for global similarity across cell types
){
  
  # 1. Run the quick pipeline for the real pathway -- does not use any global metric for transcriptome
  #    with the specified k (found previously as optimal)
  res_list = quickPipeline(pathway_genes = pathway_genes,
                           seurat_obj = seurat_obj,
                           k_final = k_final, 
                           min_genes_on = min_genes_on, 
                           min_expr = min_expr
  )
  
  # 2. Compute diversity for the cell types containing each of the pathway profiles
  # Here we use an embedding to compute the global metric (PCA, VAE, ICA, etc)
  stats_list = globalClustering(df_devel = res_list$data_frame,
                                seurat_obj = seurat_obj,
                                k_final = k_final,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric,
                                pathway_genes = pathway_genes)
  
  
  # 3. Negative control: random sets of genes from null list
  # Note: here we want two types of filtering:
  #       current: filter cell types based on the expression of the random set
  #       probably better: filter cell types that express the real pathway so the diversity is comparable
  # new option to filter the same cell types as the real pathway
  if(filter_manual){
    manual_cell_types = res_list$data_frame$cell_id
  }else{
    manual_cell_types = c()
  }
  
  # 4. Positive control: re-distribute pathway profiles in randomly selected cell types.
  # returns diversity scores only for the fake pathway (no global calculation)
  pos_control = positiveControl(df_devel=res_list$data_frame,
                                seurat_obj = seurat_obj,
                                n_random = n_samples,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric
  )
  # Merge all in data.frame
  pos_control = do.call(rbind, pos_control)
  
  
  df_diversity = rbind(
    data.frame(d = stats_list$pathway$diversity, 
               type = 'pathway'), # actual pathway
    data.frame(d = stats_list$global$diversity, 
               type = 'transcriptome'), # global dendrogram using the same cells than the real pathway
    data.frame(d = pos_control$diversity, 
               type = 'pos control')) #randomized pathway profiles across cell types
  
  df_recurrence = rbind(
    data.frame(d = stats_list$pathway$diversity * stats_list$pathway$n, 
               type = 'pathway'),
    data.frame(d = stats_list$global$diversity * stats_list$global$n, 
               type = 'transcriptome'),
    data.frame(d = pos_control$diversity * pos_control$n, 
               type ='pos control'))
  
  # We can also return the main data.frame, with class_label, diversity and rank for each profile
  clustered_data <- res_list$data_frame %>% 
    left_join(stats_list$pathway %>%
                rename(class_label = label) %>%
                dplyr::select(rank, 
                              diversity, 
                              class_label, 
                              n), 
              by = 'class_label')
  
  return(list(diversity = df_diversity, 
              recurrence = df_recurrence, 
              profiles = clustered_data, 
              rank = stats_list)
  )
}

# Computes the dispersion on cell types expressing pathway using distance in PCA space
control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                  k_final = optimal_k_pathway,
                                  seurat_obj = master_seurat,
                                  n_samples = 100,
                                  filter_manual = T,
                                  min_genes_on = min_genes_pathway, 
                                  min_expr = min_expr_threshold, 
                                  n_pcs = 100,
                                  manual_embedding = pca_proj, # PCA embedding for dispersion 
                                  dist_metric = "euclidean"
)

div_res_max_pairwise <- diverseFilt(control_res = control_res,
                                    pathway_genes = pathway_genes,
                                    diverse_quantile = diverse_quantile,
                                    type="motif"
)$diverse_mat %>% rownames()

## Median pairwise distances
globalClustering <- function(df_devel = data.frame(),   # data.frame in wide format with genes and class labels
                             pathway_genes = c(),       # List of pathway genes
                             seurat_obj = c(),          # Seurat object
                             k_final = 25,              # k value for clusters
                             n_pcs = 100,               # Number of PCs to use
                             manual_embedding = c(),    # Matrix should represent a latent space (i.e. PCA, UMAP), do not provide if using n_pcs
                             dist_metric = 'euclidean'  # Distance metric
){
  
  df_devel$class_label <- df_devel$class_label %>% as.character()
  
  # Rename df_devel by cell_id
  row.names(df_devel) <- df_devel$cell_id
  
  # Set class_label colors
  colors_1206$class_label = makeQualitativePal(length(df_devel$class_label %>% unique))
  names(colors_1206$class_label) <- df_devel$class_label %>% unique
  
  # Get the Embedding
  # umap_coords by default uses cell_id as the row.name
  if(length(manual_embedding) == 0){
    # If user specifies n_pcs
    umap_coords <- Embeddings(seurat_obj, reduction = 'pca')
    umap_coords <- umap_coords[, 1:n_pcs]
    
  }else{
    # Else use the user's Embedding, when provided directly, n_pcs is ignored
    # row.names equal cell_id
    umap_coords = manual_embedding
    
  }
  
  # Use the cells that only are in our data.frame object, with the pathway "ON"
  scaled_data = umap_coords[df_devel$cell_id, ]
  
  # Compute distance in embedding space between cells
  if(dist_metric =='euclidean')  user_dist = dist(scaled_data)  else  user_dist = dist.cosine(scaled_data)
  
  # Global clustering of cell types uses euclidean distance + ward.D2
  p_global = pheatmap(scaled_data %>% t ,
                      annotation_col = df_devel %>% 
                        dplyr::select(age, dataset, Tissue, Cell_class, class_label),
                      annotation_colors = colors_1206,
                      clustering_method = 'ward.D2',
                      cutree_cols = k_final,
                      show_rownames = F,
                      show_colnames = F,
                      fontsize = 12,
                      silent = T,
                      clustering_distance_cols = user_dist)
  
  # For each pathway label, we calculate the mean distance in the embedding space for the cell types in that label
  
  k_final = length(df_devel$class_label %>% unique)
  
  # Distance metric to be used as matrix
  # We convert the distance to a matrix so we can compute the pair-wise similarity
  dist_mat = as.matrix(user_dist)
  
  # For each class, compute the pairwise similarity
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$class_label == x]
    
    if(length(class_cells) > 1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% median())
      
    }else{      # If there is only one member in the class, then pairwise similarity is zero
      return(0)
    }
    
  }) -> diversity_pathway
  
  # Perform the same pairwise similarity computation on global labels
  
  global_labels = cutree(p_global$tree_col, 
                         k = k_final)
  df_devel$global_label = global_labels
  
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$global_label == x]
    
    if(length(class_cells)>1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% median())
      
    }else{
      return(0)
    }
    
  }) -> diversity_global
  
  # Transcriptome diversity
  df_devel %>% 
    group_by(global_label) %>% 
    count %>% 
    as.data.frame() -> global_stats
  global_stats$diversity = diversity_global
  
  # Pathway diversity
  df_devel %>% 
    group_by(class_label) %>% 
    count %>% 
    as.data.frame() %>%
    mutate(class_label = as.numeric(class_label)) %>%
    arrange(class_label) -> path_stats
  path_stats$diversity = diversity_pathway
  
  # Compute silhouette scores for this clustering
  # By default, we perform gene clustering using the cosine distance and ward.D2
  # so we also use cosine distance for the silhouette score
  
  pathway_silh = clusteringSilhouette(df = df_devel, 
                                      dist_metric = 'cosine', 
                                      pathway_genes = pathway_genes)
  
  # Join with the main data.frame
  path_stats %>% left_join(pathway_silh, by = 'class_label') -> path_stats
  # Dummy column to match later
  global_stats$mean_silh = 0
  
  path_stats$type = "pathway"
  global_stats$type = "transcriptome"
  
  path_stats %>% rename(label = class_label) -> path_stats
  global_stats %>% rename(label = global_label) -> global_stats
  
  # Rank the profiles by diversity
  path_stats$rank = rank(path_stats$diversity)
  
  return(list('pathway' = path_stats, 'global' = global_stats))
}
fullControlPathway <- function(pathway_genes = c(),      # List of pathway genes
                               k_final = c(),            # k value for clusters
                               seurat_obj = c(),         # Seurat object
                               n_samples = 100,          # No. of samples to estimate positive control
                               filter_manual = F,
                               min_expr = 0.2,           # Min. expression cutoff for gene to be "ON"
                               min_genes_on = 1,         # Min. number of genes for pathway to be "ON"
                               n_pcs = 100,              # Number of PCs to use
                               manual_embedding = c(),   # Number of principal components to consider for diversity metric
                               dist_metric = 'euclidean' # distance metric to be use in the global embedding (e.g., PCA) for global similarity across cell types
){
  
  # 1. Run the quick pipeline for the real pathway -- does not use any global metric for transcriptome
  #    with the specified k (found previously as optimal)
  res_list = quickPipeline(pathway_genes = pathway_genes,
                           seurat_obj = seurat_obj,
                           k_final = k_final, 
                           min_genes_on = min_genes_on, 
                           min_expr = min_expr
  )
  
  # 2. Compute diversity for the cell types containing each of the pathway profiles
  # Here we use an embedding to compute the global metric (PCA, VAE, ICA, etc)
  stats_list = globalClustering(df_devel = res_list$data_frame,
                                seurat_obj = seurat_obj,
                                k_final = k_final,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric,
                                pathway_genes = pathway_genes)
  
  
  # 3. Negative control: random sets of genes from null list
  # Note: here we want two types of filtering:
  #       current: filter cell types based on the expression of the random set
  #       probably better: filter cell types that express the real pathway so the diversity is comparable
  # new option to filter the same cell types as the real pathway
  if(filter_manual){
    manual_cell_types = res_list$data_frame$cell_id
  }else{
    manual_cell_types = c()
  }
  
  # 4. Positive control: re-distribute pathway profiles in randomly selected cell types.
  # returns diversity scores only for the fake pathway (no global calculation)
  pos_control = positiveControl(df_devel=res_list$data_frame,
                                seurat_obj = seurat_obj,
                                n_random = n_samples,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric
  )
  # Merge all in data.frame
  pos_control = do.call(rbind, pos_control)
  
  
  df_diversity = rbind(
    data.frame(d = stats_list$pathway$diversity, 
               type = 'pathway'), # actual pathway
    data.frame(d = stats_list$global$diversity, 
               type = 'transcriptome'), # global dendrogram using the same cells than the real pathway
    data.frame(d = pos_control$diversity, 
               type = 'pos control')) #randomized pathway profiles across cell types
  
  df_recurrence = rbind(
    data.frame(d = stats_list$pathway$diversity * stats_list$pathway$n, 
               type = 'pathway'),
    data.frame(d = stats_list$global$diversity * stats_list$global$n, 
               type = 'transcriptome'),
    data.frame(d = pos_control$diversity * pos_control$n, 
               type ='pos control'))
  
  # We can also return the main data.frame, with class_label, diversity and rank for each profile
  clustered_data <- res_list$data_frame %>% 
    left_join(stats_list$pathway %>%
                rename(class_label = label) %>%
                dplyr::select(rank, 
                              diversity, 
                              class_label, 
                              n), 
              by = 'class_label')
  
  return(list(diversity = df_diversity, 
              recurrence = df_recurrence, 
              profiles = clustered_data, 
              rank = stats_list)
  )
}

# Computes the dispersion on cell types expressing pathway using distance in PCA space
control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                  k_final = optimal_k_pathway,
                                  seurat_obj = master_seurat,
                                  n_samples = 100,
                                  filter_manual = T,
                                  min_genes_on = min_genes_pathway, 
                                  min_expr = min_expr_threshold, 
                                  n_pcs = 100,
                                  manual_embedding = pca_proj, # PCA embedding for dispersion 
                                  dist_metric = "euclidean"
)

div_res_quant_pairwise <- diverseFilt(control_res = control_res,
                                    pathway_genes = pathway_genes,
                                    diverse_quantile = diverse_quantile,
                                    type="motif"
)$diverse_mat %>% rownames()

