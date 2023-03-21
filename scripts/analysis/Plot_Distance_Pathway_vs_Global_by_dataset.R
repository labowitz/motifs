source("./scripts/analysis/imports_new.R")

pathway_list = c("Tgf-beta family receptors", # Motifs
                 "Ephrins",
                 'SRSF Splicing Protein Family',
                 'Notch receptors, Dll ligands, and Fringe proteins', 
                 "LPS and Citrate Signaling and Inflammation",     # Private
                 "NF-kB Signaling Pathway",
                 "Ras Signaling Pathway ",
                 "Scrambled Tgf-beta family receptors" # Scrambled
                 )

pathway_df <- addGenes(pathway_df = pathway_df,
                       pathway_name = "Scrambled Tgf-beta family receptors",
                       pathway_genes = genesPathway(pathway_name = "Tgf-beta family receptors", 
                                                    pathway_df = pathway_df,
                                                    seurat_obj = master_seurat))

for (j in pathway_list[[4]]){
  
  pathway_genes = genesPathway(pathway_name = j,
                               pathway_df = pathway_df,
                               seurat_obj = master_seurat)
  
  min_genes_pathway = 2 # tgfb min. number of genes expressed
  min_expr_threshold = 0.3 # tgfb minimum expression threshold for gene to be on
  optimal_k_pathway = 30 # doesn't matter for this!
  diverse_quantile = 0.9
  
  dist_metric = "cosine"
  n_pcs = 1:20
  seurat_obj = master_seurat
  
  pipe_run <- quickPipeline(seurat_obj = master_seurat,
                            pathway_genes = pathway_genes,
                            k_final = optimal_k_pathway, 
                            min_genes_on = min_genes_pathway, 
                            min_expr = min_expr_threshold
  )
  
  if (j == "Scrambled Tgf-beta family receptors"){
    pipe_run$data_frame[,pathway_genes] <- randomizeColumns(df = pipe_run$data_frame[,pathway_genes],
                                                            pathway_genes = pathway_genes)
  }
  
  print(paste(j, "cell types expressing", dim(pipe_run$data_frame)[[1]], sep=" "))
  
  # PAHTWAY DISTANCE: cosine distance in pathway gene exp.
  pathway_dist = dist.cosine(as.matrix(pipe_run$data_frame[,pathway_genes]))
  pathway_dist_df <- data.frame(t(combn(pipe_run$data_frame$cell_id,2)), 
                                as.numeric(pathway_dist))
  colnames(pathway_dist_df) <- c("cell_id_1", "cell_id_2", "pathway_dist")
  
  # GLOBAL DISTANCE: euclidean or cosine distance in PCA space
  global_coords = Embeddings(seurat_obj, reduction='pca')
  global_coords = global_coords[pipe_run$data_frame$cell_id, n_pcs] # already set row names as cell id
  global_dist = dist(global_coords)
  global_dist_df <- data.frame(t(combn(rownames(global_coords),2)), 
                               as.numeric(global_dist))
  
  colnames(global_dist_df) <- c("cell_id_1", "cell_id_2", "global_dist")
  
  merged_df <- merge(x = pathway_dist_df, 
                     y = global_dist_df, 
                     by.x = c("cell_id_1", 
                              "cell_id_2"))
  
  alpha = 0.01
  if (j == 'Notch receptors, Dll ligands, and Fringe proteins'){
    alpha = 0.1
  }
  
  p <- ggplot(merged_df, aes(pathway_dist, global_dist)) + 
    geom_point(alpha = alpha) + 
    ggtitle(paste("Pathway vs. Global Distance: ", j, sep="")) + 
    xlab("Pathway Distance (cosine similarity on gene exp.)") + 
    ylab("Global Distance (euclidean distance on 20 PCs)")
  
  ggsave(paste("./scripts/figures/Pathway_Global_Distance/euclidean_global_",j,".pdf", sep = ""))
  ggsave(paste("./scripts/figures/Pathway_Global_Distance/euclidean_global_",j,".png", sep = ""))
  
  
  for (i in unique(pipe_run$data_frame$dataset)){
    
    pathway_dist = dist.cosine(as.matrix(pipe_run$data_frame[pipe_run$data_frame$dataset == i,pathway_genes]))
    
    pathway_dist_df <- data.frame(t(combn(pipe_run$data_frame[pipe_run$data_frame$dataset == i,]$cell_id,2)), 
                                  as.numeric(pathway_dist))
    
    colnames(pathway_dist_df) <- c("cell_id_1", "cell_id_2", "pathway_dist")
    
    global_coords = Embeddings(seurat_obj, reduction='pca')
    global_coords = global_coords[pipe_run$data_frame[pipe_run$data_frame$dataset == i,]$cell_id, n_pcs] # already set row names as cell id
    global_dist = dist(global_coords)
    global_dist_df <- data.frame(t(combn(rownames(global_coords),2)), 
                                 as.numeric(global_dist))
    
    colnames(global_dist_df) <- c("cell_id_1", "cell_id_2", "global_dist")
    
    merged_df <- merge(x = pathway_dist_df, 
                       y = global_dist_df, 
                       by.x = c("cell_id_1", 
                                "cell_id_2"))
    
    alpha = 0.5
    if (j == 'Notch receptors, Dll ligands, and Fringe proteins'){
      alpha = 1
    }
    if (i == "E6.5_8.5_Chan"){
      alpha=0.1
    }
    
    if (j == 'SRSF Splicing Protein Family'){
      alpha=0.1
    }
    p <- ggplot(merged_df, aes(pathway_dist, global_dist)) + 
      geom_point(alpha = alpha) + 
      ggtitle(paste("Pathway vs. Global Distance:", j, i, sep=" ")) + 
      xlab("Pathway Distance (cosine similarity on gene exp.)") + 
      ylab("Global Distance (euclidean distance on 20 PCs)")
    
    ggsave(paste("./scripts/figures/Pathway_Global_Distance/euclidean_global_",j,"_", i,".pdf", sep = ""))
    ggsave(paste("./scripts/figures/Pathway_Global_Distance/euclidean_global_",j,"_", i,".png", sep = ""))
    
  }
}

