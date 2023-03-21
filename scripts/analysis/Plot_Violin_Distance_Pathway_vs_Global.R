source("./scripts/analysis/imports_new.R")

pathway_name =  "Tgf-beta family receptors" # tgfb pathway name
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj = master_seurat)

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
optimal_k_pathway = 30 # optimal pathway components, computed from z-score
diverse_quantile = 0.9

dist_metric = "cosine"
n_pcs = 1:20
seurat_obj = master_seurat

## To generate a DataFrame with the gene exp. of cell states with pathway "ON"
pipe_run <- quickPipeline(seurat_obj = master_seurat,
                          pathway_genes = pathway_genes,
                          k_final = optimal_k_pathway, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold
)

# PAHTWAY DISTANCE: cosine distance in pathway gene exp.
pathway_dist = dist.cosine(as.matrix(pipe_run$data_frame[pathway_genes]))
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

global_dist_df %>% 
  ggplot(aes(x=global_dist)) +
  geom_histogram(aes(x=global_dist,y=after_stat(density))) +
  geom_density() -> p

pg <- ggplot_build(p)

#peak_indx <- c(F,diff(sign(c(diff(pg$data[[2]]$y))))==-2,F) %>% which()
#topcounts <- pg$data[[2]]$x[peak_indx]

peak_indx <- c(F,diff(sign(c(diff(pg$data[[2]]$y))))==2,F) %>% which()
botcounts <- pg$data[[2]]$x[peak_indx]

global_dist_df$global_dist <- findInterval(global_dist_df$global_dist, append(c(0), botcounts))

global_dist_df$global_dist <- append(c(0), botcounts)[match(global_dist_df$global_dist, global_dist_df$global_dist %>% unique())] %>% as.factor()

merged_df <- merge(x = pathway_dist_df, 
                   y = global_dist_df, 
                   by.x = c("cell_id_1", 
                            "cell_id_2"))

p <- p + 
  geom_vline(xintercept = append(c(0), botcounts),linetype="dashed") + 
  xlab("Global Distance (euclidean distance on 20 PCs)") + 
  theme_classic()

ggsave(paste("./scripts/figures/Pathway_Global_Distance/violin_cutoff_hist_",pathway_name,".pdf", sep = ""))
ggsave(paste("./scripts/figures/Pathway_Global_Distance/violin_cutoff_hist_",pathway_name,".png", sep = ""))

p <- ggplot(merged_df, aes(pathway_dist, global_dist)) + 
  geom_violin()+
  ggtitle("Pathway vs. Global Distance") + 
  xlab("Pathway Distance (cosine similarity on gene exp.)") + 
  ylab("Global Distance (euclidean distance on 20 PCs)") + 
  theme_classic()

ggsave(paste("./scripts/figures/Pathway_Global_Distance/violin_plot_",pathway_name,".pdf", sep = ""))
ggsave(paste("./scripts/figures/Pathway_Global_Distance/violin_plot_",pathway_name,".png", sep = ""))
