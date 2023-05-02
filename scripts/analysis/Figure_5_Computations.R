source("./scripts/analysis/imports.R")
source("./scripts/analysis/Figure_5_Functions.R") # Script with Fig. 5 functions

## Run the silhouette score and dispersion computations for all pathways
# Directories to save results to
silh_res_dir = "./scripts/figures/peak_analysis/silhouette_res/silh_rds/"
dispersion_dir = "./scripts/figures/peak_analysis/dispersion/disp_rds/"

min_genes_pathway = 2
min_expr_threshold = 0.3
diverse_quantile = 0.9

max_k = 200 # Number of clusters for which we want to compute silhoutte scores

runPipeline <- function(pathway_name =""){

  print( paste( "Running ..", pathway_name , " " , Sys.time() ))
  
  pathway_genes <- genesPathway(pathway_name = pathway_name,
                                pathway_df = pathway_df,
                                seurat_obj = master_seurat)
  
  if (length(pathway_genes) > 7){
    
    # if it founds an error, it won't break but it will save the pathway as NA
    # we can then read the outputs and re-run the NA with a lower number of maxk 
    # the most likely cause of error is too many clusters 
    silh_plt = tryCatch({silhouettePlot(pathway_genes = pathway_genes, 
                                        min_genes_on = min_genes_pathway, 
                                        min_expr = min_expr_threshold, 
                                        n_bootstraps = 100,
                                        seurat_obj = master_seurat,
                                        max_k = max_k) 
    }, error = function(e) { print(e); return ("NA")
    }, finally = { print("error handled") } )
    
    # save the file! This is what we want!! 
    if(silh_plt != "NA"){
      saveRDS(silh_plt, paste(silh_res_dir, pathway_name, ".RDS", sep=""))
      
      silh_z_plt <- silhouette_zscore(silh_result = silh_plt,
                                      min_expr = min_expr_threshold,
                                      x_offset = 6,
                                      max.y = 0.43,
                                      min.y = 0.17,
                                      k_max = silh_plt[[2]]$k %>% max()
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
      
      # make the plot -- less important 
      print( paste( "Done ..", pathway_name , " " , Sys.time() ))
    }else{
      print( paste( "Error ..", pathway_name , " " , Sys.time() ))
    }
    
  }

}

mclapply(pathway_df$pathway %>% unique(), runPipeline, mc.cores = detectCores()-1)

# For pathways that failed, try to compute silhouette for fewer clusters
max_k = 150 # Number of clusters for which we want to compute silhoutte scores

runPipeline <- function(pathway_name =""){
  
  if (!file.exists(paste(silh_res_dir, pathway_name, ".RDS", sep=""))){
    
    print( paste( "Running ..", pathway_name , " " , Sys.time() ))
    
    pathway_genes <- genesPathway(pathway_name = pathway_name,
                                  pathway_df = pathway_df,
                                  seurat_obj = master_seurat)
    
    if (length(pathway_genes) > 7){
      
      # if it founds an error, it won't break but it will save the pathway as NA
      # we can then read the outputs and re-run the NA with a lower number of maxk 
      # the most likely cause of error is too many clusters 
      silh_plt = tryCatch({silhouettePlot(pathway_genes = pathway_genes, 
                                          min_genes_on = min_genes_pathway, 
                                          min_expr = min_expr_threshold, 
                                          n_bootstraps = 100,
                                          seurat_obj = master_seurat,
                                          max_k = max_k) 
      }, error = function(e) { print(e); return ("NA")
      }, finally = { print("error handled") } )
      
      # save the file! This is what we want!! 
      if(silh_plt != "NA"){
        saveRDS(silh_plt, paste(silh_res_dir, pathway_name, ".RDS", sep=""))
        
        silh_z_plt <- silhouette_zscore(silh_result = silh_plt,
                                        min_expr = min_expr_threshold,
                                        x_offset = 6,
                                        max.y = 0.43,
                                        min.y = 0.17,
                                        k_max = silh_plt[[2]]$k %>% max()
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
        
        # make the plot -- less important 
        print( paste( "Done ..", pathway_name , " " , Sys.time() ))
      }else{
        print( paste( "Error ..", pathway_name , " " , Sys.time() ))
      }
      
    }
    
  }
  
}

mclapply(pathway_df$pathway %>% unique(), runPipeline, mc.cores = detectCores()-1)
