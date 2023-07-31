source("./scripts/analysis/imports.R")
output_dir = "./scripts/analysis/outputs/aws_pathways/"

pathway_list = c("CXCR4 Signaling Pathway",
                 "Lysophosphatidic Acid LPA6 Signalling",
                 "Rac 1 Cell Motility Signaling Pathway",
                 "Insulin Signalling",
                 "Cadmium Induces DNA Synthesis and Proliferation in Macrophages ",
                 "Apoptotic DNA Fragmentation and Tissue Homeostasis",
                 "GnRH Signaling Pathway",
                 "Growth Hormone Signaling Pathway ",
                 "Ubiquitin–Proteasome Pathway",
                 "RNA-splicing by SR protein family",
                 "Eph-Ephrin"
)


perc_k_finder <- function(z_score = c(),
                          percentile = 0.7
){
  idx_max_z = which.max(z_score)
  max_z = max(z_score, na.rm=T)
  
  vals = rev(z_score[idx_max_z:length(z_score)])
  
  ret_idx = which.min(abs(vals-(max_z * percentile)))
  
  return(length(vals) - ret_idx + idx_max_z)
  
}

## Modified z-score plot function that overestimates the k_max
## Computes and plots Z-score based on bootstrapped pathway and randomized expression values
silhouette_zscore <- function(silh_result,                        # Result from running silhouettePlot function
                              min_expression = min_expr_threshold,# Min. expression cutoff for gene to be "ON"
                              x_offset = 1,                       # X-axis of k (number of clusters) begins plotting at x_offset + 1
                              pathway_name = 'Notch',             # Pathway name
                              min.y = 0.1,                        # Min. cut-off on y-axis of the silhouette score plot
                              max.y = 0.7                         # Max. cut-off on y-axis of the silhouette score plot
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
  
  df = data.frame(k=1:100, z = z_score)
  df %>% dplyr::filter(k>x_offset) -> df
  
  # Plots
  ## Plot the pathway genes and randomized genes bootstrapped values
  p1 <- boot_df %>% ggplot(aes(x = k , y = m)) +
    geom_ribbon(aes(ymin = m - s, ymax = m + s) ,fill = 'gray', alpha = 0.2) +
    geom_line(color = 'black', size = 1) + geom_ribbon(data = boot_df_control, aes(ymin = m - s, ymax = m + s), alpha = 0.2, fill = 'gray') +
    geom_line(data = boot_df_control, color = 'gray', size = 1)  + theme_pubr(base_size = 10) + ylab('Silhouette') +
    ggtitle(pathway_name) + scale_y_continuous(limits = c(min.y, max.y))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank())
  
  ## Plot the Z-score
  p2 <- ggplot(data = df, aes(x = k, y = z)) + geom_line(data = df, aes(x = k, y = z), color = rgb(54, 55, 149, max = 255), size = 1) +
    geom_vline(xintercept = perc_k_finder(z_score), linetype = "dashed") + theme_pubr(base_size = 10) +
    ylab("Z-score") + xlab("Number clusters") + theme(axis.title.y = element_text(colour = , color = rgb(54, 55, 149, max = 255)))
  
  # Merge the two plots into 1 using ggarrange
  g <- ggarrange(p1, NULL, p2, ncol = 1, align = "hv", heights = c(1, -0.175, 0.5))
  
  return(g)
}

make_plts <- function(pathway_name){
  silh_result = paste(output_dir, pathway_name, "_silh_plt.RDS", sep = "")
  
  silhouette_zscore(silh_result = readRDS(silh_result),
                    min_expression = min_expr_threshold,
                    pathway_name = pathway_name, 
                    min.y = -0.05, 
                    max.y = 0.55
  )
}

grob_list = lapply(pathway_list, make_plts)
g = arrangeGrob(grobs = grob_list)
ggsave(paste(output_dir, "silh_plts.pdf", sep = ""), g, height = 36, width = 16)