source("./scripts/analysis/imports.R")

# Function to compute the cosine distance of two gene expression profiles
computePathwayDist <- function(pipe_run){
  
  pathway_dist = dist.cosine(as.matrix(pipe_run$data_frame[,pathway_genes]))
  
  pathway_dist_df <- data.frame(t(combn(pipe_run$data_frame$cell_id,2)), 
                                as.numeric(pathway_dist))
  
  colnames(pathway_dist_df) <- c("cell_id_1", "cell_id_2", "pathway_dist")
  
  return(pathway_dist_df)
  
}

# Computes the pathway, global, and null dist. of profile distances in expressing cell types
plotPathwayGlobalDist <- function(pathway_genes = c(),
                                  pathway_df = pathway_df,
                                  seurat_obj = master_seurat,
                                  k_final = optimal_k_pathway,
                                  min_genes_on = min_genes_pathway,
                                  min_expr = min_expr_threshold){
  
  pipe_run <- quickPipeline(seurat_obj = seurat_obj,
                            pathway_genes = pathway_genes,
                            k_final = optimal_k_pathway, 
                            min_genes_on = min_genes_on, 
                            min_expr = min_expr
  )
  
  scrambled_pipe_run <- pipe_run
  
  scrambled_pipe_run$data_frame[,pathway_genes] <- randomizeColumns(df = pipe_run$data_frame[,pathway_genes],
                                                                    pathway_genes = pathway_genes)
  
  real_pathway_dist_df <- computePathwayDist(pipe_run)
  
  scrambled_pathway_dist_df <- computePathwayDist(scrambled_pipe_run)
  colnames(scrambled_pathway_dist_df)[colnames(scrambled_pathway_dist_df) == 'pathway_dist'] <- 'random_dist'

  # GLOBAL DISTANCE: euclidean or cosine distance in PCA space
  global_coords = Embeddings(seurat_obj, reduction='pca')
  global_coords = global_coords[pipe_run$data_frame$cell_id, n_pcs] # already set row names as cell id
  
  global_dist = dist(global_coords)
  
  global_dist_df <- data.frame(t(combn(rownames(global_coords),2)), 
                               as.numeric(global_dist))
  
  colnames(global_dist_df) <- c("cell_id_1", "cell_id_2", "global_dist")
  
  df <- merge(x = real_pathway_dist_df, 
               y = scrambled_pathway_dist_df,
               by.x = c("cell_id_1", 
                        "cell_id_2"))
  
  df <- merge(x = df, 
              y = global_dist_df, 
              by.x = c("cell_id_1", 
                       "cell_id_2"))
  
  return(df)
}

pathway_name = "Tgf-beta family receptors"

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
optimal_k_pathway = 30 # doesn't matter for this!
diverse_quantile = 0.9

n_pcs = 1:30

pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj = master_seurat)

plotPathwayGlobalDist(pathway_genes = pathway_genes) -> df

# Generates a list of dataframes with 100 scrambled estimates for null dist.
df_list <- lapply(rep(list(pathway_genes), 100), plotPathwayGlobalDist)
