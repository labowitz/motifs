source("./scripts/analysis/imports.R")
source("./scripts/analysis/Figure_5_Functions.R") # Script with Fig. 5 functions

## Run the silhouette score and dispersion computations for all pathways
# Directories to save results to
silh_res_dir = "./scripts/figures/peak_analysis/silhouette_res/silh_rds/"
dispersion_dir = "./scripts/figures/peak_analysis/dispersion/disp_rds/"

min_genes_pathway = 2
min_expr_threshold = 0.3
diverse_quantile = 0.9

for (pathway_name in unique(pathway_df$pathway)){
  
  pathway_genes <- genesPathway(pathway_name = pathway_name,
                                pathway_df = pathway_df,
                                seurat_obj = master_seurat)
  
  if (length(pathway_genes) >= 7){
    silh_plt = silhouettePlot(pathway_genes = pathway_genes, 
                              min_genes_on = min_genes_pathway, 
                              min_expr = min_expr_threshold, 
                              n_bootstraps = 100,
                              seurat_obj = master_seurat)
    
    saveRDS(silh_plt, paste(silh_res_dir, pathway_name, ".RDS", sep=""))
    
    silh_z_plt <- silhouette_zscore(silh_result = silh_plt,
                                    min_expr = min_expr_threshold,
                                    x_offset = 6,
                                    max.y = 0.43,
                                    min.y = 0.17
    )
    
    optimal_k_pathway <- perc_k_finder(z_score = silh_z_plt[[2]],
                                       percentile = 0.9)
    
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
    
    saveRDS(control_res, paste(dispersion_dir, pathway_name, ".RDS", sep=""))
  }
  
}

