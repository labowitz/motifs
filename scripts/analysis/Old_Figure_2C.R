source("./scripts/analysis/imports_new.R")
output_dir <- "./scripts/figures/"

pathways = c("Bmp_Tgfb", "Notch", "Eph_r", "Wnt")
min_gene_list = c(2,2,2,2)
min_expr_list = c(0.2, 0.2, 0.3, 0.3)
optimal_k_list = c(30, 30, 54, 30)

for (i in 1:length(pathways)){
  
  pathway_df = readRDS("./data/processed_data/all_pathways.RDS")
  
  pathway_name =  pathways[[i]]
  min_genes_pathway = min_gene_list[[i]]
  min_expr_threshold = min_expr_list[[i]]
  optimal_k_pathway = optimal_k_list[[i]]
  
  pathway_genes = genesPathway(pathway_name = pathway_name,
                               pathway_df = pathway_df)
  
  control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                    k_final = optimal_k_pathway,
                                    seurat_obj = master_seurat, # seurat object
                                    null_list = hvg_genes, #list of highly variable genes 
                                    n_samples = 100, 
                                    filter_manual = T,
                                    min_genes_on = min_genes_pathway, 
                                    min_expr = min_expr_threshold, 
                                    n_pcs = 100, # how many PCs to use
                                    manual_embedding = pca_proj, # PCA embedding for diversity 
                                    dist_metric = "euclidean"
  )
  
  g_umap <- global_umap(control_res = control_res,
                        seurat_obj = master_seurat,
                        use_pca = T,
                        n_pcs = 1:20,
                        clust_method = "ward.D2",
                        dist_metric = "cosine")
  
  ## Add the pathway statistics
  g_umap <- g_umap+ggtitle()
  
  ggsave(paste(output_dir, "Figure_2C_", i, ".pdf", sep=""))
  
}