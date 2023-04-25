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

# Optimal number of pathway components (computed previously with silh. score)
optimal_k_pathway = 30

# Make a plot of clustered profiles for cell types with pathway "ON"
pdf(paste(fig_dir, "Figure_3B.pdf",sep=""))
pipe_run <- quickPipeline(seurat_obj = master_seurat,
                          pathway_genes = pathway_genes,
                          k_final = optimal_k_pathway, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold
)
print(pipe_run$plot)
dev.off()

  # Computing the optimal silhouette score
silh_plt = silhouettePlot(pathway_genes = pathway_genes, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold, 
                          n_bootstraps = 100,
                          seurat_obj = master_seurat)

silh_z_plt <- silhouette_zscore(silh_result = silh_plt,
                                min_expr = min_expr,
                                x_offset = 6,
                                max.y = 0.43,
                                min.y = 0.17
)
silh_z_plt
ggsave(paste(fig_dir, "Figure_3A.pdf", sep = ""), width = 7, height = 7)

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

# Plot global dendrogram with pathway states
pdf(paste(fig_dir, "Figure_3C.pdf", sep=""))
glob_dendr <- global_dendr(control_res = control_res,
                           seurat_obj = master_seurat,
                           use_pca = T,
                           n_pcs = 1:20,
                           clust_method = 'ward.D2',
                           dist_metric ='cosine')
print(glob_dendr$plt)
dev.off()

# Ranking of dispersion scores of profiles
rank_plot <- rank_diversity(pathway_genes = pathway_genes,
                            min_genes_on = min_genes_pathway,
                            dist_metric = "euclidean",
                            make_plot = T,
                            k_final = optimal_k_pathway, 
                            min_expr = min_expr_threshold, 
                            manual_embedding = pca_proj,
                            seurat_obj = master_seurat
)
ggsave(paste(fig_dir, "Figure_4B.pdf", sep = ""), width = 7, height = 7)

# ECDF of dispersion scores
ecdf_plot = ecdf_diversity(control_res)
ggsave(paste(fig_dir, "Figure_4C.pdf", sep = ""), width = 7, height = 7)

# Motif heatmap, motifs defined as profiles above 90th percentile of dispersion
pdf(paste(fig_dir, "Figure_4D.pdf", sep=""))
motif_heat <- motif_heatmap(control_res = control_res,
                            pathway_genes = pathway_genes,
                            diverse_quantile = diverse_quantile,
                            type="motif"
)
print(motif_heat)
dev.off()

# Motif distribution across tissues and organs.
motif_ct <- motif_ct_heatmap(control_res = control_res,
                             pathway_genes = pathway_genes,
                             diverse_quantile = diverse_quantile,
                             type="motif"
)
ggsave(plot = motif_ct, paste(fig_dir, "Figure_4E.pdf", sep=""))

# Private profiles heatmap
pdf(paste(fig_dir, "Figure_S3D.pdf", sep=""))
priv_heat <- motif_heatmap(control_res = control_res,
                            pathway_genes = pathway_genes,
                            diverse_quantile = diverse_quantile,
                            type="private"
)
print(priv_heat)
dev.off()

# UMAP showing cell types with pathway "ON"
g_umap <- global_umap(control_res = control_res,
                      seurat_obj = master_seurat,
                      use_pca = T,
                      n_pcs = 1:20,
                      clust_method = "ward.D2",
                      dist_metric = "cosine")

ggsave(plot = g_umap, paste(fig_dir, "Figure_2C.pdf", sep=""))

# Bar graph of cell types with gene "ON."
counts_bar <- geneCounts(seurat_obj = master_seurat,
                         pathway_genes = pathway_genes,
                         min_genes_on = min_genes_pathway,
                         min_expr = min_expr_threshold)
ggsave(paste(fig_dir, "Figure_S2A.pdf", sep=""))

# ECDF of number of genes "ON" at different thresholds.
ecdf_thresh <- ecdfThresh(seurat_obj = master_seurat, 
                          pathway_genes = pathway_genes,
                          min_genes_on = min_expr_threshold)

ggsave(plot = ecdf_thresh, filename = paste(fig_dir, 
                                            "Figure_S2B.pdf", 
                                            sep=""))

# Heatmap of no. of profiles co-expressing genes pairwise above threshold.
coexp_heat <- coexpHeatmap(seurat_obj = master_seurat,
                           pathway_genes = pathway_genes,
                           min_genes_on = min_genes_pathway,
                           min_expr = min_expr_threshold)

ggsave(plot = coexp_heat, filename = paste(fig_dir, 
                                           "Figure_S2C.pdf", 
                                           sep=""))
