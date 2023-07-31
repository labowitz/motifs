## Directories for the files
source("./scripts/analysis/imports.R")
output_dir = "./scripts/analysis/outputs/aws_pathways/"

# Import data
data_real_plot <- read.table("./data/raw_data/allPathways_listGenes_dec2021.tsv", header = T, sep = "\t")

# Correct pathway names
data_real_plot$pathway[data_real_plot$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
data_real_plot$pathway[data_real_plot$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
data_real_plot$pathway[data_real_plot$pathway=='Srsf'] <- 'RNA-splicing by SR protein family'
data_real_plot$pathway[data_real_plot$pathway=='Eph_r'] <- 'Eph A-B receptors'
data_real_plot$pathway[data_real_plot$pathway=='Wnt'] <- 'Frizzled and Lrp5/6 receptors for Wnt B/Catenin Signaling'
data_real_plot$pathway[data_real_plot$pathway=='Fgfr'] <- 'FGF cell signaling proteins'

# Revised genesPathway function to grab genes frmo data_real_plot object instead
genesPathway <- function(which_pathway = 'Notch'){
  
  this_pathway = data_real_plot %>% dplyr::filter(pathway == which_pathway) %>% pull(gene)
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]
  
  return(this_pathway)
}

# Revised quickPipeline function to grab genes frmo data_real_plot object instead
quickPipeline <- function(master_seurat = c(),              # Seurat object
                          which_pathway = 'Notch',          # Pathway name to draw from all_pathways object
                          k_final = 25,                     # k value for clusters
                          min_genes_on = 1,                 # Min. number of genes for pathway to be "ON"
                          min_expr = 0.2,                   # Min. expression cutoff for gene to be "ON"
                          which_profiles = 'devel',         # Include developmental, adult, or both datasets
                          rand_control = F,                 #
                          silent_plot = T,                  # Whether or not to return heatmap of data and class labels
                          manual_filter_cells = c(),        #
                          verbose = F,
                          save_pdf = F,
                          pdf_file = ""
){
  if(!rand_control){
    this_pathway = data_real_plot %>% dplyr::filter(pathway == which_pathway) %>% pull(gene)
  }else{
    this_pathway = which_pathway # random set directly as input
  }
  
  # Get pathway genes in the dataset
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]
  
  # Get the MinMax scaled and normalized data and annotations
  df_devel <- normalizedDevel(master_seurat = master_seurat,
                              this_pathway,
                              sat_val = 0.99,
                              fill_zero_rows = F,
                              which_datasets = which_profiles)
  
  # Add a column of cell_type based on the different column values
  df_devel %>% mutate(cell_type = paste(global_cluster, '--',Tissue,': ', cell_ontology_class,'-', age, sep="")) -> df_devel
  
  # Compute the number of genes on in the pathway
  df_devel$genes_on = rowSums(df_devel[, this_pathway] > min_expr)
  
  row.names(df_devel) <- df_devel$cell_type
  
  # Two types of filtering.
  # Either by minimum expression of pathway genes OR user-specified list of cell types.
  if(length(manual_filter_cells) == 0){
    # Filter out cells with pathway "ON"
    df_devel %>% dplyr::filter(genes_on > min_genes_on) -> df_devel
  }else{
    # Else filter user-specified cell types
    df_devel %>% dplyr::filter(cell_id %in% manual_filter_cells) -> df_devel
    
    # minority of gene sets will have non-expressing cells
    # this happens very rarely and only 1 or 2 cells will have no expression, so we can remove them
    expressing_cells = df_devel$cell_id[rowSums(df_devel[, this_pathway]) > 0]
    df_devel %>% dplyr::filter(cell_id %in% expressing_cells) -> df_devel
  }
  
  # Heatmap to for computing the distance tree for cell types
  p = pheatmap(df_devel[,this_pathway ],
               #annotation_row = df_devel %>% select(dataset, Cell_class),
               #annotation_colors = colors_1206,
               show_rownames = T,
               fontsize = 5,
               cutree_rows = k_final,
               clustering_method = 'ward.D2',
               clustering_distance_rows = dist.cosine(as.matrix(df_devel[, this_pathway])),
               cluster_cols = F,
               silent = T)
  
  # Get class labels
  cos_labels = cutree(p$tree_row, k = k_final) # Cosine tree at specified k-cut
  df_devel$class_label = cos_labels # Map class labels
  
  df_devel$class_label <- as.character(df_devel$class_label)
  
  # Map colors to the class_labels and store in list called annoCols
  annoCols <- makeQualitativePal(k_final, glasbey_use = T, skip = 1) # Make Glasbey palette with k_final categorical colors
  names(annoCols) <- unique(df_devel$class_label)
  annoCols <- list("class_label" = annoCols)
  
  # Heatmap to return with class_label annotations
  if(!save_pdf){
    p = pheatmap(df_devel[, this_pathway],
                 show_rownames = F,
                 fontsize = 5,
                 annotation_row = df_devel %>% select(class_label),
                 annotation_colors = annoCols,col = black_pal(100),
                 cutree_rows = k_final, clustering_method = 'ward.D2',
                 clustering_distance_rows = dist.cosine(as.matrix(df_devel[, this_pathway])),
                 cluster_cols = F,
                 silent = silent_plot)
  }else{
    p = pheatmap(df_devel[, this_pathway],
                 show_rownames = F,
                 fontsize = 5,
                 annotation_row = df_devel %>% select(class_label),
                 annotation_colors = annoCols,
                 col = black_pal(100),
                 cutree_rows = k_final,
                 clustering_method = 'ward.D2',
                 clustering_distance_rows = dist.cosine(as.matrix(df_devel[, this_pathway])),
                 cluster_cols = F,
                 silent = silent_plot,
                 filename = pdf_file,
                 height = 5,
                 width = 5)
  }
  
  df_devel$class_label <- as.numeric(df_devel$class_label)
  
  df_devel %>% gather(key = 'gene', value ='expression', this_pathway) %>%
    group_by(class_label, gene) %>%
    summarise(mean_expr = mean(expression), n_cell_types = n()) %>%
    spread(gene, mean_expr) %>% tibble::column_to_rownames(var = "class_label") -> x
  
  # Return count matrix, annotations, and heatmap for selected pathway
  return(list('matrix' = x, 'df_devel' = df_devel, 'profiles' = p))
}

# List of pathways to do silhouette analysis on

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
                 "Eph A-B receptors"
                 )

# Parameters to show
min_expr_threshold <- 0.3
min_genes_pathway <- 2

for(i in pathway_list){
  pathway_name <- i
  
  silh_plt = silhouettePlot(which_pathway = pathway_name, 
                            min_ON = min_genes_pathway, 
                            min_expr = min_expr_threshold, 
                            n_bootstraps = 100)
  
  saveRDS(silh_plt, paste(output_dir, i, "_silh_plt.RDS", sep = ""))
  
  g <- silhouette_zscore(silh_result = silh_plt,
                         pathway_name = pathway_name,
                         x_offset = 6,
                         max.y = 0.6,
                         min.y = 0.1 # adjust axis parameters 
  )
  g
  ggsave(paste(output_dir, i, "_silh_z-score.pdf", sep = ""), width = 7, height = 7)
}

perc_k_finder <- function(z_score = c(),
                          percentile = 0.7
){
  idx_max_z = which.max(z_score)
  max_z = max(z_score, na.rm=T)
  
  vals = rev(z_score[idx_max_z:length(z_score)])
  
  ret_idx = which.min(abs(vals-(max_z * percentile)))
  
  return(length(vals) - ret_idx + idx_max_z)
  
}

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
    geom_vline(xintercept = which.max(z_score), linetype = "dashed") + theme_pubr(base_size = 10) +
    ylab("Z-score") + xlab("Number clusters") + theme(axis.title.y = element_text(colour = , color = rgb(54, 55, 149, max = 255)))
  
  # Merge the two plots into 1 using ggarrange
  g <- ggarrange(p1, NULL, p2, ncol = 1, align = "hv", heights = c(1, -0.175, 0.5))
  
  return(list(z_score, g))
}

df_kvals = data.frame(name = pathway_list, k = rep(0, length(pathway_list)))

for (i in 1:length(pathway_list)){
  silh_result = paste(output_dir, pathway_list[[i]], "_silh_plt.RDS", sep = "")
  
  res = silhouette_zscore(silh_result = readRDS(silh_result),
                          min_expression = min_expr_threshold,
                          pathway_name = pathway_list[[i]], 
                          min.y = -0.05, 
                          max.y = 0.55
  )
  
  ggsave(paste(output_dir, pathway_list[[i]], "_silh_z-score.pdf", sep = ""), width = 7, height = 7)
  
  df_kvals$k[[i]] = perc_k_finder(z_score = res[[1]])                    
}