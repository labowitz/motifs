source("./scripts/analysis/imports_new.R")

# Run pipeline in parallel
runSilhBoot <- function(pathway_name = "",
                        output_dir = "",
                        pathway_df = pathway_df,
                        min_genes_on = 2,
                        min_expr = 0.3,
                        n_bootstraps = 5,
                        max_k = 100,
                        seurat_obj = c()){
  
  print( paste( "Running ..", pathway_name , " " , Sys.time() ))
  
  pathway_genes = genesPathway(pathway_name = pathway_name,
                               pathway_df = pathway_df)
  
  # if it founds an error, it won't break but it will save the pathway as NA
  # we can then read the outputs and re-run the NA with a lower number of maxk 
  # the most likely cause of error is too many clusters 
  silh_plt = tryCatch({silhouettePlot(pathway_genes = pathway_genes, 
                                      min_genes_on = min_genes_on,
                                      min_expr = min_expr, 
                                      n_bootstraps = n_bootstraps,
                                      seurat_obj = seurat_obj,
                                      max_k = max_k)
    
  }, error = function(e) { print(e); return ("NA")
  }, finally = { print("error handled") } )
  
  # save the file! This is what we want!! 
  if(silh_plt != "NA"){
    
    saveRDS(silh_plt, paste(output_dir, 
                            pathway_name, 
                            "_silh_plt.RDS", 
                            sep = ""))
  }else{
    print( paste( "Error ..", pathway_name , " " , Sys.time()))
  }
  
}

min_genes_on = 2
min_expr = 0.3
seurat_obj = master_seurat
max_k = 200
output_dir = "data/processed_data/Silhouette_PathBank/"

mclapply(unique(pathway_df$pathway), # List of pathways
         FUN = runSilhBoot, 
         min_genes_on=min_genes_on,
         pathway_df,
         min_expr=min_expr,
         seurat_obj = seurat_obj,
         max_k = max_k,
         output_dir = output_dir,
         mc.cores = 8) 
