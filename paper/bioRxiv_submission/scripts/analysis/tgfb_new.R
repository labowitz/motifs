library(ggplotify)
source("./scripts/analysis/imports_trial.R")

## Change these variables
filename = "tgfb_trial" # how name is appended to files for pathway

## Directories for the files
output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")
fig_dir = "./scripts/figures/tgfb_trial/"

pathway_df = readRDS("./data/processed_data/all_pathways.RDS")

pathway_name =  "Bmp_Tgfb" # tgfb pathway name
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df)

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
optimal_k_pathway = 30 # optimal pathway components, computed from z-score
diverse_quantile = 0.9

## Make a plot of k_opt profiles for cell types with pathway "ON
pipe_run <- quickPipeline(seurat_obj = master_seurat,
                          pathway_genes = pathway_genes,
                          k_final = optimal_k_pathway, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold
)

## Computing the optimal silhouette score
silh_plt = silhouettePlot(pathway_genes = pathway_genes, 
                          min_genes_on = min_genes_on, 
                          min_expr = min_expr, 
                          n_bootstraps = 5,
                          seurat_obj = master_seurat)

silh_z_plt <- silhouette_zscore(silh_result = silh_plt,
                                min_expr = min_expr,
                                x_offset = 6,
                                max.y = 0.43,
                                min.y = 0.17
)
silh_z_plt
ggsave(paste(fig_dir, "silh_new.pdf", sep = ""), width = 7, height = 7)

## Running fullControlPathway
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

## Plotting
# ECDF plots
ecdf_plot = ecdf_diversity(control_res)
ggsave(paste(fig_dir, "ecdf_div_new.pdf", sep = ""), width = 7, height = 7)

# Diversity rank plots
rank_plot <- rank_diversity(pathway_genes = pathway_genes,
                            min_genes_on = min_genes_pathway,
                            dist_metric = "euclidean",
                            make_plot = T,
                            k_final = optimal_k_pathway, 
                            min_expr = min_expr_threshold, 
                            manual_embedding = pca_proj,
                            seurat_obj = master_seurat
)
ggsave(paste(fig_dir, "rank_div_new.pdf", sep = ""), width = 7, height = 7)

# Plot global dendrogram with pathway states

pdf(paste(fig_dir, "glob_dendr_new.pdf", sep=""))
glob_dendr <- global_dendr(control_res = control_res,
                           seurat_obj = master_seurat,
                           use_pca = T,
                           n_pcs = 1:20,
                           clust_method = 'ward.D2',
                           dist_metric ='cosine')
dev.off()

# 6. Tissue and organ distributions of selected profiles
pdf(paste(fig_dir, "motif_heat_new.pdf", sep=""))
motif_heat <- motif_heatmap(control_res = control_res,
                            pathway_genes = pathway_genes,
                            diverse_quantile = diverse_quantile,
                            type="motif"
)
dev.off()

motif_ct <- motif_ct_heatmap(control_res = control_res,
                             pathway_genes = pathway_genes,
                             diverse_quantile = diverse_quantile,
                             type="motif"
)
ggsave(plot = motif_ct, paste(fig_dir, "motif_ct_new.pdf", sep=""))

g_umap <- global_umap(control_res = control_res,
                      seurat_obj = master_seurat,
                      use_pca = T,
                      n_pcs = 1:20,
                      clust_method = "ward.D2",
                      dist_metric = "cosine")

ggsave(plot = g_umap, paste(fig_dir, "global_umap_new.pdf", sep=""))

p_umap <- profiles_umap(control_res = control_res,
                        seurat_obj = master_seurat,
                        use_pca = T,
                        n_pcs = 1:20,
                        clust_method = "ward.D2",
                        dist_metric = "cosine",
                        ncol = 6,
                        text_size = 10)

ggsave(plot = p_umap, filename = paste(fig_dir, "profile_umap_new.pdf", sep=""))

counts_bar <- geneCounts(seurat_obj = master_seurat,
                         pathway_genes = pathway_genes,
                         min_genes_on = min_genes_pathway,
                         min_expr = min_expr_threshold)
ggsave(paste(fig_dir, "counts_barplot_new.pdf", sep=""))

coexp_heat <- coexpHeatmap(seurat_obj = master_seurat,
                           pathway_genes = pathway_genes,
                           min_genes_on = min_genes_pathway,
                           min_expr = min_expr_threshold)

ggsave(plot = coexp_heat, filename = paste(fig_dir, 
                                           "coexp_heatmap_new.pdf", 
                                           sep=""))
split_pathway_genes = split(pathway_genes,
                            cut(seq_along(pathway_genes), 
                                3, 
                                labels = FALSE))

grobs = list(textGrob(paste("\n",
                            paste(split_pathway_genes$`1`, 
                                  collapse="\n"),
                            sep=""),
                      just = "center",
                      gp=gpar(fontsize=14)),
             textGrob(paste("Pathway genes\n",
                            paste(split_pathway_genes$`2`, 
                                  collapse="\n"),
                            sep=""), 
                      just = "center",
                      gp=gpar(fontsize=14)),
             textGrob(paste("\n",
                            paste(split_pathway_genes$`3`, 
                                  collapse="\n"),
                            sep=""), 
                      just = "center",
                      gp=gpar(fontsize=14)), 
             coexp_heat +
               theme(plot.margin = unit(c(0,0,0,2), "cm")), 
             as.ggplot(pipe_run$plot), 
             silh_z_plt[[1]] + 
               theme(axis.title = element_text(size=16), 
                     plot.title = element_text(size=20),
                     axis.text = element_text(size=14)),
             #counts_bar + 
             #theme(axis.title = element_text(size=20), 
             #plot.title = element_text(size=20)), 
             g_umap + 
               theme(plot.title = element_text(size=20),
                     plot.margin = unit(c(0,1,1.5,0), "cm")), 
             ecdf_plot + 
               theme(axis.title = element_text(size=16), 
                     plot.title = element_text(size=20),
                     axis.text = element_text(size=14),
                     legend.text = element_text(size=10),
                     legend.title = element_text(size=12),
                     plot.margin = unit(c(0,1,0,0), "cm")), 
             rank_plot + 
               theme(axis.title = element_text(size=16), 
                     plot.title = element_text(size=20),
                     axis.text = element_text(size=14),
                     plot.margin = unit(c(0,1,0,0), "cm")),
             as.ggplot(motif_heat$plot) + 
               ggtitle(label=paste("Motifs, Dispersion >= ", 100*diverse_quantile, "th percentile", sep="")) + 
               theme(plot.title = element_text(size=20)),
             motif_ct + 
               theme(plot.margin = unit(c(2,0,0,0), "cm")) + 
               ggtitle("Motif Tissue Composition") +
               theme(plot.title = element_text(size=20))
)

lay = rbind(c(1,2,3,6,6,6,5,5),
            c(4,4,4,6,6,6,5,5),
            c(4,4,4,6,6,6,5,5),
            #c(6,6,6,NA,NA,NA,5,5),
            #c(6,6,6,NA,NA,NA,5,5),
            #c(7,7,8,8,9,9,5,5),
            c(7,7,8,8,9,9,5,5),
            c(7,7,8,8,9,9,5,5),
            c(10,10,10,11,11,11,5,5),
            c(10,10,10,11,11,11,5,5),
            c(10,10,10,11,11,11,5,5),
            c(10,10,10,NA,NA,NA,5,5))

m.grid.arrange <- function(p1, 
                           lay1, 
                           p2,
                           k_final,
                           pathway_name = c()
) 
{
  pdf(file = NULL) #invisible
  
  ml <- list(grid.arrange(grobs = p1, 
                          layout_matrix = lay1,
                          vp=viewport(width=0.975, height=0.975) # adds space between page and plot
  ), 
  p2
  )
  
  return(marrangeGrob(grobs=ml,
                      nrow=1,
                      ncol=1,
                      top=textGrob(paste(pathway_name,
                                         " (k_opt = ", optimal_k_pathway, ")",
                                         sep=""), 
                                   gp=gpar(fontsize=36)))
  )
  dev.off() #invisible
}

ml <- m.grid.arrange(p1 = grobs, 
                     lay1 = lay,
                     p2 = p_umap,
                     k_final = optimal_k_pathway,
                     pathway_name = pathway_name
)

ggsave(filename = paste(fig_dir, "trial_arrange.pdf", sep=""), plot = ml, width = 17, height = 22)
