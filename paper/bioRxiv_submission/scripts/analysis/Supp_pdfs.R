source("./scripts/analysis/imports_new.R")

silh_res_dir = "scripts/figures/peak_analysis/silhouette_res/silh_rds/"

# Import data -- lists of genes 
pathway_df <- read.table("./data/raw_data/allPathways_listGenes_dec2021.tsv", 
                         header = T, 
                         sep = "\t")
# Correct pathway names
pathway_df$pathway[pathway_df$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
pathway_df$pathway[pathway_df$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
pathway_df$pathway[pathway_df$pathway=='Srsf'] <- 'RNA-splicing by SR protein family'
pathway_df$pathway[pathway_df$pathway=='Eph_r'] <- 'Eph A-B receptors'
pathway_df$pathway[pathway_df$pathway=='Wnt'] <- 'Frizzled and Lrp5 6 receptors for Wnt B Catenin Signaling'
pathway_df$pathway[pathway_df$pathway=='Fgfr'] <- 'FGF cell signaling proteins'
# Add a combined version of Eph ligands + receptors 
pathway_df <- rbind(pathway_df, data.frame(pathway = 'Eph receptors and ligands', 
                                           gene = pathway_df %>% 
                                             dplyr::filter(grepl(pattern = 'Eph', pathway)) %>% 
                                             dplyr::pull(gene))) 
# fix the names of some pathways 
pathway_df$pathway <- pathway_df$pathway %>% str_replace('/', ' ')
pathway_df$pathway <- pathway_df$pathway %>% str_replace('\\(', ' ')
pathway_df$pathway <- pathway_df$pathway %>% str_replace('\\)', ' ')

fig_5_df <- read.csv("scripts/figures/peak_analysis/silhouette_res/dispersion_states_figure5.csv")

pathway_df <- pathway_df[pathway_df$pathway %in% fig_5_df$pathway_name,]

# Functions:
perc_k_finder <- function(z_score = c(),
                          percentile = 0.9
){
  idx_max_z = which.max(z_score)
  max_z = max(z_score, na.rm=T)
  
  vals = rev(z_score[idx_max_z:length(z_score)])
  
  ret_idx = which.min(abs(vals-(max_z * percentile)))
  
  return(length(vals) - ret_idx + idx_max_z)
  
}

silhouette_zscore <- function(silh_result,                        # Result from running silhouettePlot function
                              min_expr = 0.2,                     # Min. expression cutoff for gene to be "ON"
                              x_offset = 1,                       # X-axis of k (number of clusters) begins plotting at x_offset + 1
                              min.y = 0.1,                        # Min. cut-off on y-axis of the silhouette score plot
                              max.y = 0.7,                         # Max. cut-off on y-axis of the silhouette score plot
                              k_max = 100,
                              percentile_peak = 0.9
){
  # Obtain results from running the silhouettePlot function
  # Format of silh_result is a list with silh_result = (plot, boot_df, boot_df_control)
  x = silh_result
  
  # Extract data.frames from list and remove offset on the x-axis
  boot_df = x[[2]] %>% dplyr::filter(k>x_offset)
  boot_df_control = x[[3]] %>% dplyr::filter(k>x_offset)
  
  # Compute the Z-score, which is the difference between the pathway and randomized silhouettes divided by the randomized std. dev.
  diff = x[[2]]$m-x[[3]]$m;
  z_score = diff/x[[3]]$s
  max_z = max(z_score, na.rm = T)
  
  df = data.frame(k=1:k_max, z = z_score)
  df %>% dplyr::filter(k>x_offset) -> df
  
  # Plots
  ## Plot the pathway genes and randomized genes bootstrapped values
  p1 <- boot_df %>% ggplot(aes(x = k , y = m)) +
    geom_ribbon(aes(ymin = m - s, ymax = m + s), 
                fill = 'gray', 
                alpha = 0.2) +
    geom_line(color = 'black', 
              linewidth = 1) + 
    geom_ribbon(data = boot_df_control, 
                aes(ymin = m - s, 
                    ymax = m + s), 
                alpha = 0.2, 
                fill = 'gray') +
    geom_line(data = boot_df_control, 
              color = 'gray', 
              linewidth = 1) + 
    theme_pubr(base_size = 16) + 
    ylab('Silhouette') + 
    ggtitle("Silhouette and Z-score") + 
    scale_y_continuous(limits = c(min.y, 
                                  max.y)) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.line.x = element_blank()
    )
  
  ## Plot the Z-score
  p2 <- ggplot(data = df, 
               aes(x = k, 
                   y = z)) + 
    geom_line(data = df, 
              aes(x = k, 
                  y = z), 
              color = rgb(54, 55, 149, max = 255), 
              linewidth = 1) +
    geom_vline(xintercept = perc_k_finder(z_score, 
                                          percentile_peak), 
               linetype = "dashed") + 
    theme_pubr(base_size = 16) +
    ylab("Z-score") + 
    xlab("Number clusters") + 
    theme(axis.title.y = element_text(color = rgb(54, 55, 149, max = 255)))
  
  # Merge the two plots into 1 using ggarrange
  g <- ggarrange(p1, 
                 NULL, 
                 p2, 
                 ncol = 1, 
                 align = "hv", 
                 heights = c(1, -0.175, 0.5)
  )
  
  return(list(g, z_score))
}

## Save pathway data.frame
saved_files <- paste(silh_res_dir, 
                     list.files(path=silh_res_dir,
                                pattern = "_ale.RDS"), 
                     sep ="") # ale is my personal ID 

# if the pathway finished and the file was saved correctly: 
full_pathway_list <- pathway_df %>% 
  group_by(pathway) %>% 
  count %>% 
  as.data.frame %>% 
  dplyr::filter(n>7) %>% 
  pull(pathway)

saved_idx <- lapply(full_pathway_list, 
                    FUN = function(x) sum(grepl(x, saved_files))) %>% 
  unlist %>% 
  as.logical

# load only the pathways that actually finished 
pathway_list_dispersion <- full_pathway_list[saved_idx]

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

generatePlots(seurat_obj = master_seurat,
              pathway_name = pathway_list_dispersion[[i]],
              pathway_df = pathway_df,
              min_genes_pathway = 2,
              min_expr_threshold = 0.3,
              diverse_quantile = 0.9,
              percentile_peak = 0.9,
              silh_res_dir = silh_res_dir)

generatePlots <- function(seurat_obj = c(),
                          pathway_name = "",
                          pathway_df = data.frame(),
                          min_genes_pathway = 2, # eph min. number of genes expressed
                          min_expr_threshold = 0.3, # eph minimum expression threshold for gene to be on
                          diverse_quantile = 0.9,
                          percentile_peak = 0.9,
                          silh_res_dir = ""
                          ){
  
  pathway_genes = genesPathway(pathway_name = pathway_name,
                               pathway_df = pathway_df,
                               seurat_obj = seurat_obj)
  
  ## Computing the optimal silhouette score
  silh_result = readRDS(paste(silh_res_dir,
                              pathway_name,
                              "_silh_plt_ale.RDS",
                              sep = ""))
  
  silh_z_plt <- silhouette_zscore(silh_result = silh_result,
                                  min_expr = min_expr,
                                  x_offset = 6,
                                  max.y = 0.43,
                                  min.y = 0.17,
                                  k_max = silh_result[[2]]$k %>% 
                                    max()
  )
  
  optimal_k_pathway = perc_k_finder(silh_z_plt[[2]], 
                                    percentile_peak)
  
  pipe_run <- quickPipeline(seurat_obj = master_seurat,
                            pathway_genes = pathway_genes,
                            k_final = optimal_k_pathway, 
                            min_genes_on = min_genes_pathway, 
                            min_expr = min_expr_threshold
  )
  
  ## Running fullControlPathway
  control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                    k_final = optimal_k_pathway,
                                    seurat_obj = seurat_obj, # seurat object
                                    null_list = hvg_genes, #list of highly variable genes 
                                    n_samples = 100, 
                                    filter_manual = F,
                                    min_genes_on = min_genes_pathway, 
                                    min_expr = min_expr_threshold, 
                                    n_pcs = 100, # how many PCs to use
                                    manual_embedding = pca_proj, # PCA embedding for diversity 
                                    dist_metric = "euclidean"
  )
  
  ecdf_plot = ecdf_diversity(control_res)
  
  rank_plot <- rank_diversity(pathway_genes = pathway_genes,
                              min_genes_on = min_genes_pathway,
                              dist_metric = "euclidean",
                              make_plot = T,
                              k_final = optimal_k_pathway, 
                              min_expr = min_expr_threshold, 
                              manual_embedding = pca_proj,
                              seurat_obj = master_seurat
  )
  
  motif_heat <- motif_heatmap(control_res = control_res,
                              pathway_genes = pathway_genes,
                              diverse_quantile = diverse_quantile,
                              type="motif"
  )
  
  motif_ct <- motif_ct_heatmap(control_res = control_res,
                               pathway_genes = pathway_genes,
                               diverse_quantile = diverse_quantile,
                               type="motif"
  )
  
  g_umap <- global_umap(control_res = control_res,
                        seurat_obj = master_seurat,
                        use_pca = T,
                        n_pcs = 1:20,
                        clust_method = "ward.D2",
                        dist_metric = "cosine")
  
  p_umap <- profiles_umap(control_res = control_res,
                          seurat_obj = master_seurat,
                          use_pca = T,
                          n_pcs = 1:20,
                          clust_method = "ward.D2",
                          dist_metric = "cosine",
                          ncol = 6,
                          text_size = 10)
  
  coexp_heat <- coexpHeatmap(seurat_obj = master_seurat,
                             pathway_genes = pathway_genes,
                             min_genes_on = min_genes_pathway,
                             min_expr = min_expr_threshold)
  
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
              c(7,7,8,8,9,9,5,5),
              c(7,7,8,8,9,9,5,5),
              c(10,10,10,11,11,11,5,5),
              c(10,10,10,11,11,11,5,5),
              c(10,10,10,11,11,11,5,5),
              c(10,10,10,NA,NA,NA,5,5))
  
  ml <- m.grid.arrange(p1 = grobs, 
                       lay1 = lay,
                       p2 = p_umap,
                       k_final = optimal_k_pathway,
                       pathway_name = pathway_name
  )
  
  ggsave(filename = paste(silh_res_dir, pathway_name, ".pdf", sep=""), plot = ml, width = 17, height = 22)
  
}

