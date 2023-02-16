## Directories for the files
source("./imports.R")
output_dir = "./"


## "C:/Users/bacte/Documents/Caltech/rnaseq/Niv/fig5/min_exp_update/output"
#setwd("C:/Users/bacte/Documents/Caltech/rnaseq/Niv/fig5/min_exp_update/")
output_dir_ = 'output/'


# Import data -- lists of genes 
data_real_plot <- read.table("./data/raw_data/allPathways_listGenes_dec2021.tsv", header = T, sep = "\t")

# Correct pathway names
data_real_plot$pathway[data_real_plot$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
data_real_plot$pathway[data_real_plot$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
data_real_plot$pathway[data_real_plot$pathway=='Srsf'] <- 'RNA-splicing by SR protein family'
data_real_plot$pathway[data_real_plot$pathway=='Eph_r'] <- 'Eph A-B receptors'
data_real_plot$pathway[data_real_plot$pathway=='Wnt'] <- 'Frizzled and Lrp5 6 receptors for Wnt B Catenin Signaling'
data_real_plot$pathway[data_real_plot$pathway=='Fgfr'] <- 'FGF cell signaling proteins'

# Add a combined version of Eph ligands + receptors 
data_real_plot <- rbind(data_real_plot, data.frame(pathway = 'Eph receptors and ligands', gene = data_real_plot %>% dplyr::filter(grepl(pattern = 'Eph', pathway )) %>% pull(gene ) ) ) 

# fix the names of some pathways 
data_real_plot$pathway  <- data_real_plot$pathway %>% str_replace('/', ' ')
data_real_plot$pathway  <- data_real_plot$pathway %>% str_replace('\\(', ' ')
data_real_plot$pathway  <- data_real_plot$pathway %>% str_replace('\\)', ' ')

# Revised genesPathway function to grab genes frmo data_real_plot object instead
genesPathway <- function(which_pathway = 'Notch'){
  
  this_pathway = data_real_plot %>% dplyr::filter(pathway == which_pathway) %>% pull(gene)
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]
  
  return(this_pathway)
}

# Revised quickPipeline function to grab genes from data_real_plot object instead
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
# Make list of pathways 
data_real_plot %>% group_by(pathway) %>% count() %>% dplyr::filter(n>7) %>% pull(pathway) -> filter_pathway_list 


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
                 "Eph A-B receptors",
                 'Frizzled and Lrp5 6 receptors for Wnt B Catenin Signaling',
                 'Notch receptors, Dll ligands and Fringe proteins',
                 'Tgf-beta family receptors', 
                 'Eph receptors and ligands'
                 )

# # Parameters to show
min_expr_threshold <- 0.3
min_genes_pathway <- 2

# for(i in pathway_list){
#   pathway_name <- i
  
#   silh_plt = silhouettePlot(which_pathway = pathway_name, 
#                             min_ON = min_genes_pathway, 
#                             min_expr = min_expr_threshold, 
#                             n_bootstraps = 10)
  
#   saveRDS(silh_plt, paste(output_dir_, i, "_silh_plt.RDS", sep = ""))
  
#   g <- silhouette_zscore(silh_result = silh_plt,
#                          pathway_name = pathway_name,
#                          x_offset = 6,
#                          max.y = 0.6,
#                          min.y = 0.1 # adjust axis parameters 
#   )
#   g
#   ggsave(paste(output_dir, i, "_silh_z-score.pdf", sep = ""), width = 7, height = 7)
# }


### END silh RUN 
################
###############
###############

# Run pipeline in parallel
runPipeline <- function(pathway_name ="" ){
    
    print( paste( "Running ..", pathway_name , " " , Sys.time() ))
    
		# if it founds an error, it won't break but it will save the pathway as NA
		# we can then read the outputs and re-run the NA with a lower number of maxk 
		# the most likely cause of error is too many clusters 
		silh_plt = tryCatch({silhouettePlot(which_pathway = pathway_name, 
                          min_ON = min_genes_pathway, 
                          min_expr = min_expr_threshold, max_k = 200, 
                          n_bootstraps = 10)  
        }, error = function(e) { print(e); return ("NA")
        }, finally = { print("error handled") } )
    
    # save the file! This is what we want!! 
    if(silh_plt != "NA"){
	    saveRDS(silh_plt, paste(output_dir, pathway_name, "_silh_plt_ale.RDS", sep = ""))
    # make the plot -- less important 
			print( paste( "Done ..", pathway_name , " " , Sys.time() ))
		}else{
   	  print( paste( "Error ..", pathway_name , " " , Sys.time() ))
		}
    
}

## Optimal K 
## Dispersion calculations


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
                              min_expression = min_expr_threshold,# Min. expression cutoff for gene to be "ON"
                              x_offset = 1,                       # X-axis of k (number of clusters) begins plotting at x_offset + 1
                              pathway_name = 'Notch',             # Pathway name
                              min.y = 0.1,                        # Min. cut-off on y-axis of the silhouette score plot
                              max.y = 0.7,                        # Max. cut-off on y-axis of the silhouette score plot
                              percentile_peak =0.7
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
  
  df = data.frame(k=1:length(z_score), z = z_score)
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
    geom_vline(xintercept = which.max(z_score), linetype = "dashed") + 
    geom_vline(xintercept = perc_k_finder(z_score,percentile_peak ), linetype = "dashed", color = "red") + 
    
    theme_pubr(base_size = 10) +
    ylab("Z-score") + xlab("Number clusters") + theme(axis.title.y = element_text(colour = , color = rgb(54, 55, 149, max = 255)))
  
  # Merge the two plots into 1 using ggarrange
  g <- ggarrange(p1, NULL, p2, ncol = 1, align = "hv", heights = c(1, -0.175, 0.5))
  
  return(list(z_score, g))
}





#######################
# Niv Update Nov 5 2022
# 0. Run the silhouette functions for all pathways and save results ... DONE 
# 1. Read RDS output
# 2. Find the optimal number of clusters  --- calls perc_k_finder() 
# 3. Run the pipeline and compute dispersion 
process_silh <- function(save_file = F){
  # returns a data.frame of pathways and optimal number of clusters 
  output_dir_ = "./"
  saved_files <- paste(output_dir_, list.files()[grepl('ale.RDS', list.files() )], sep ="") # ale is my personal ID 

  # if the pathway finished and the file was saved correctly: 
  full_pathway_list <- data_real_plot %>% group_by(pathway) %>% count %>% as.data.frame %>% dplyr::filter(n>7) %>% pull(pathway)
  saved_idx <- lapply(full_pathway_list, FUN = function(x) sum(grepl(x, saved_files )) )  %>% unlist %>% as.logical

  # load only the pathways that actually finished 
  pathway_list_dispersion <- full_pathway_list[saved_idx]

  # prepare the data frame for pathways with output 
  df_kvals = data.frame(name = pathway_list_dispersion, k = rep(0, length(pathway_list_dispersion)))
  output_dir_ = "./"
  # Read output from all pathwways and run the pipeline using 
  # the optimal number of clusters found above 
  for (i in 1:length(pathway_list_dispersion)){
    silh_result = paste(output_dir_, pathway_list_dispersion[[i]], "_silh_plt_ale.RDS", sep = "")
    
    res = silhouette_zscore(silh_result = readRDS(silh_result),
                            min_expression = min_expr_threshold,
                            pathway_name = pathway_list_dispersion[[i]], 
                            min.y = -0.05, 
                            max.y = 0.55
    )
    
    #ggsave(paste(output_dir, pathway_list[[i]], "_silh_z-score_retrial.pdf", sep = ""), width = 7, height = 7)
    
    df_kvals$k[[i]] = perc_k_finder(z_score = res[[1]])                    
  }

  if(save_file)
    saveRDS(df_kvals, paste(output_dir, "df_kvals.RDS"))


  return(df_kvals)
}

# Run the pipeline <- let's do this part in parallel 
### 



# 1. For each pathway, run the pipeline using the optimal number of clusters identified by process_silh() 
# 2. Export a data.frame of pathway, k , dispersion 
# 3. Save the full output as .RDS for each pathway 
output_dir_dispersion ="dispersion/"
computeDipersions <- function(
                              use_min_k = F, # whether to choose the optimal k by considering the minimum value of k that is close to the peak (or the max -- default)
                              overwrite_files = F #  if we are running only new pathways and we don't want to re-run pathways that finished already 
){

  # Find all the pathway results 
  saved_files <- paste(output_dir_, list.files()[grepl('ale.RDS', list.files() )], sep ="") # ale is my personal ID 

  # if the pathway finished and the file was saved correctly: 
  full_pathway_list <- data_real_plot %>% group_by(pathway) %>% count %>% as.data.frame %>% dplyr::filter(n>7) %>% pull(pathway)
  saved_idx <- lapply(full_pathway_list, FUN = function(x) sum(grepl(x, saved_files )) )  %>% unlist %>% as.logical

  # load only the pathways that actually finished (silhoutte)
  pathway_list_dispersion <- full_pathway_list[saved_idx]


  # subset only pathways for which we don't have DISPERSION results already 
  if(!overwrite_files){

    #check which pathways have dispersion results already (if we don't want to overwrite data)
    done_pathways <- list.files(output_dir_dispersion)
    # exclude them from the list of pahtways to run 
    done_idx <- sapply(pathway_list_dispersion, function(x){ sum(grepl(x, done_pathways)) }) %>% unlist %>% as.logical
    pathway_list_dispersion <- pathway_list_dispersion[!done_idx]
  }


  for (i in 1:length(pathway_list_dispersion)){

    if(!use_min_k){
      optimal_k_pathway = df_kvals %>% dplyr::filter(name == pathway_list_dispersion[i])  %>% pull(k)
    }else{
      optimal_k_pathway = df_kvals %>% dplyr::filter(name == pathway_list_dispersion[i])  %>% pull(min_k_peak)
    }


    print( paste( "Running ", pathway_list_dispersion[i], " k= ", optimal_k_pathway, Sys.time() ))
    control_res  = fullControlPathway(this_pathway = pathway_list_dispersion[[i]],
                                    k_pathway = optimal_k_pathway , 
                                    filter_pathway = 'both', # adult + devel 
                                    this_seurat = master_seurat, # seurat object
                                    null_list = hvg_genes, #list of highly variable genes 
                                    n_samples = 100, 
                                    filter_manual = T,
                                    min_expr_gene = min_expr_threshold, 
                                    min_genes_ON = min_genes_pathway, # default > 2 genes ON 
                                    n_pcs_global = 100, # how many PCs to use
                                    embedding_matrix = pca_proj # PCA embedding for diversity 
                                    )
    
    saveRDS(control_res, paste(output_dir_dispersion, pathway_list_dispersion[[i]], "_diversity.RDS", sep=""))

    print(paste("DONE ", pathway_list_dispersion[[i]], Sys.time() ))
  }


}


# New code to compute dispersion 
parseDispersion <- function(pathway_list_dispersion = list() 
                            ){
  dispersion_list = list() 

  for(p in pathway_list_dispersion){
      
      
      path_disp = readRDS(paste0( 'dispersion/', p, '_diversity.RDS'))
      
      df_path = path_disp[[4]]$pathway
      df_path$pathway_name  = p
      dispersion_list[[p]] = df_path
  }

  dispersion_stats = do.call(rbind, dispersion_list)

  return(dispersion_stats )
}
# make statistics 
plot_dispersion_distribution <- function(
                      dispersion_stats = data.frame() , 
                      labels_idx = c(9, 10, 12,13,14),
                      min_silh = 0.3,                   # threshold for considering only clusters with high silhouette score 
                      min_n_cluster = 2,                # do not include clusters with only one cell type 
                      df_kvals = data.frame(),           # data.frame with optimal number of clusters 
                      use_mean = F,
                      include_metabolic = F, 
                      use_palette  = 'Set2'
){
  #min_silh = 0.1 # only consider clusters with positive silhouette score -- otherwise we could bias the dispersion calculations to bad clusters 
  df_plot <- dispersion_stats %>% 
    dplyr::filter(mean_silh > min_silh, n> min_n_cluster) %>% 
    group_by(pathway_name) %>% 
    summarise( 
              silh = mean(mean_silh), 
              median_disp = median(diversity) , n_clusters = n() )


  select_pathways = pathway_list[labels_idx ]


  df_plot_labeled <- df_plot %>% left_join(df_kvals %>% rename(pathway_name = name ), by ="pathway_name") %>%
    mutate( pathway_label  = ifelse(pathway_name %in% select_pathways,str_trunc(pathway_name, 12) ,"")) %>% 
    dplyr::filter(k > 3)

  df_plot_labeled <- df_plot_labeled %>% 
                        mutate(pathway_class = ifelse(grepl("bolis|ycle|Gluconeogenesis|ynthesis|Ammonia|Succinate|Fructose|Lysine|Pentose|Fatty Acids", pathway_name),
                                         "metabolism","signaling/protein") )

  # LPA pathway appears 3 times 
  df_plot_labeled <- df_plot_labeled %>% dplyr::filter(!grepl(pattern = "LPA4|LPA5|LPA1|LPA3|LPA6", pathway_name)) 


  # merge with number of gnees 
  df_plot_labeled %>% left_join(data_real_plot %>% rename(pathway_name = pathway) %>% 
                      group_by(pathway_name) %>% summarise(n_genes = n() ) , by ="pathway_name") %>% 
                      mutate(k_gene_ratio = k/n_genes) -> df_plot_labeled

  if(!include_metabolic)
    df_plot_labeled <- df_plot_labeled %>% dplyr::filter(pathway_class != 'metabolism')

 
  g <- df_plot_labeled %>%
      ggplot(aes(x = median_disp, y = k_gene_ratio, size = silh, label = pathway_label, color = pathway_class)) + 
                          geom_point() + geom_text()

  return(list(g, df_plot_labeled ) )
}


plotly_dispersion_distribution <- function(
                      dispersion_stats = data.frame() , 
                      labels_idx = c(9, 10, 12,13,14),
                      min_silh = 0.1,                   # threshold for considering only clusters with high silhouette score 
                      min_n_cluster = 1,                # do not include clusters with only one cell type 
                      df_kvals = data.frame(),           # data.frame with optimal number of clusters 
                      use_mean = T
){
  #min_silh = 0.1 # only consider clusters with positive silhouette score -- otherwise we could bias the dispersion calculations to bad clusters 
  df_plot <- dispersion_stats %>% 
    dplyr::filter(mean_silh > min_silh, n> min_n_cluster) %>% 
    group_by(pathway_name) %>% 
    summarise(
              silh = mean(mean_silh), 
              median_disp = median(diversity) , n_clusters = n() )


  select_pathways = pathway_list[labels_idx ]


  df_plot_labeled <- df_plot %>% left_join(df_kvals %>% rename(pathway_name = name ), by ="pathway_name") %>%
    mutate( pathway_label  = ifelse(pathway_name %in% select_pathways,str_trunc(pathway_name, 8) ,"")) %>% 
    dplyr::filter(k > 3)

  df_plot_labeled <- df_plot_labeled %>% 
                        mutate(pathway_class = ifelse(grepl("bolis|ycle|Gluconeogenesis|Fructose|Lysine|Pentose|Fatty Acids", pathway_name),
                                         "metabolism","signaling/protein") )

  # LPA pathway appears 3 times 
  df_plot_labeled <- df_plot_labeled %>% dplyr::filter(!grepl(pattern = "LPA2|LPA4", pathway_name)) 


  # merge with number of gnees 
  df_plot_labeled %>% left_join(data_real_plot %>% rename(pathway_name = pathway) %>% 
                      group_by(pathway_name) %>% summarise(n_genes = n() ) , by ="pathway_name") %>% 
                      mutate(k_gene_ratio = k/n_genes) -> df_plot_labeled

  fig <- plot_ly(data = df_plot_labeled, x = ~median_disp, y = ~k_gene_ratio, color = ~pathway_class, 
              colors = 'Set1',
              text = df_plot_labeled$pathway_name,
              hoverinfo = 'text',
              
              marker = list(size = 10,
                             #color = 'rgba(255, 182, 193, .9)',
                             line = list(#color = 'rgba(152, 0, 0, .8)',
                                         width = 2)))

  return(list(fig, df_plot_labeled ) )
}



#res = plot_dispersion_distribution(dispersion_stats , labels_idx = c(9,10,12,13,14,1,3), min_silh = 0.3, min_n_cluster = 2,  df_kvals = df_kvals , use_mean = F, include_metabolic = F, use_palette = 'Paired')
 
#res[[1]] + theme(text = element_text(size = 20)) + ylab('N clusters / n genes') + xlab('Dispersion')


# Compute peak width scores for all pathways 

# Peak score function 
plot_peakdness_score <- function(z_score_raw  = c() , 
                                 pathway_name = "", 
                                 fil_win = 5, 
                                 percentile = 0.9,  # % of max to use as the peak height -- smaller values result in broader peaks 
                                 peak_width_threshold  = 0.1){
    par(mfrow = c(2,1))
    
    # smooth the z-score 
    z_score <- stats::filter(z_score_raw, rep(1,fil_win)/fil_win, sides = 1) 
    # absolute max 
    idx_max_z = which.max(z_score)
    max_z = max(z_score, na.rm=T)
    
    
    plot(z_score_raw, type = "l", main = pathway_name )
    
    lines(z_score , type = "l", col = "red")
    abline(v = idx_max_z)
    abline(h = max_z * percentile, col ="blue")
    # plot the optimal k as detected by our percentile function 
    abline(v = perc_k_finder(z_score_raw, percentile) ,lty = 2, col = "blue")  # only for visualization
    
    #subplot 2 -- relative difference to the peak percentile 
    plot(abs(  z_score - (max_z * percentile))/max_z,  type = "l")
    abline(h = peak_width_threshold  , col ="red")
    
    # absolute distance of z-score values relative to the peak 
    peak_dist <- abs(  z_score - (max_z * percentile))/max_z
    
    # After applying the smoothing filter, we get some NA values 
    peak_dist_vals <- peak_dist[!is.na(peak_dist)]
    # how many k values have a z-score within 10% of the peak 
    # if too many k values have a z-score close to the max (peak) then we can classify 
    # this peak as "broad" 
    peak_width <- sum(peak_dist_vals <= peak_width_threshold )
    peak_width <- peak_width/length(peak_dist_vals ) # this should account for pathways with k_max < 200 (default for most pathways)
    
    # let's save the range values 
    # min k value with high silhouette 
    min_k <- which(peak_dist_vals <= peak_width_threshold)[1] + fil_win
    # max value in the range 
    max_k <- rev(which(peak_dist_vals <= peak_width_threshold) )[1] + fil_win  
    
    # let's add vertical lines for the range of the peak 
    abline(v= min_k, col = 'blue')
    abline(v = max_k, col = 'blue') 
    return(list(peak_width, min_k, max_k) ) # retuns a value after plotting 
}

# Will generate and save plots of peaks for all pathways 
# Returns a data.frame with the width scores and the optimal k for all pathways 
# NOTE: pathways with low peak scores should not be considered as their "optimal k" is not a good approximation 
peak_width_scores <- function(
                              use_percentile = 0.95, # threshold for considering values within the peak range (% of max z-score)
                              smooth_window = 3,     # filtering window size for smoothing the z-score data before finding the max value. 
                              pct_opt = 0.9,         # will be depracted!!: : percentile of max z-score that will be used to identify the optimal k - this can be different from the above percentile 
                              peak_threshold = 0.1   # relative distance from the max to consider a k value within the peak 
                              ){
  output_dir_ = "./"
  saved_files <- paste(output_dir_, list.files()[grepl('ale.RDS', list.files() )], sep ="") # ale is my personal ID 

  # if the pathway finished and the file was saved correctly: 
  full_pathway_list <- data_real_plot %>% group_by(pathway) %>% count %>% as.data.frame %>% dplyr::filter(n>7) %>% pull(pathway)
  saved_idx <- lapply(full_pathway_list, FUN = function(x) sum(grepl(x, saved_files )) )  %>% unlist %>% as.logical

  # load only the pathways that actually finished 
  pathway_list_dispersion <- full_pathway_list[saved_idx]

  # prepare the data frame for pathways with output 
  df_kvals = data.frame(name = pathway_list_dispersion, 
                        k = rep(0, length(pathway_list_dispersion)), 
                        width = rep(0, length(pathway_list_dispersion)))
  output_dir_ = "./"
  # Read output from all pathwways and run the pipeline using 
  # the optimal number of clusters found above 
  for (i in 1:length(pathway_list_dispersion)){
      silh_result = paste(output_dir_, pathway_list_dispersion[[i]], "_silh_plt_ale.RDS", sep = "")
      
      res = silhouette_zscore(silh_result = readRDS(silh_result),
                              min_expression = min_expr_threshold,
                              pathway_name = pathway_list_dispersion[[i]], 
                              min.y = -0.05, 
                              max.y = 0.55
      )
      
      pdf(paste('./peak_analysis/', pathway_list_dispersion[[i]], '_silhPeak.pdf', sep = ""))
      
      # Computes the peak width 
      peak_width_list  <- plot_peakdness_score(z_score_raw = res[[1]], 
                                        pathway_name = pathway_list_dispersion[[i]], 
                                        fil_win = smooth_window , percentile = use_percentile, peak_width_threshold = peak_threshold)
      
      dev.off() 
      df_kvals$k[[i]] = perc_k_finder(z_score = res[[1]], percentile = pct_opt)  # 0.9 is the default value --- used in the manuscript 
      df_kvals$width[[i]] <- peak_width_list[[1]]
      df_kvals$min_k_peak[i] <- peak_width_list[[2]]
      df_kvals$max_k_peak[i] <- peak_width_list[[3]]
      
  }


  return(df_kvals) 


} # WE can now merge df_kvals with the dispersion scores 


# PIPELINE #############################################################
########################################################################
# run the above function 
# 1. optimal K 
# df_kvals <- peak_width_scores(use_percentile = 0.95, smooth_window = 3, pct_opt = 0.9 , peak_threshold = 0.1)
# # optional: if using min_k::::df_kvals <- df_kvals %>% mutate(k_midrange = k, k = min_k_peak ) 
# 
# Plot: Visualize the width distribution 
# df_kvals$width %>% hist(breaks =20 , col = 'lightgray', border = 'black', xlab = 'Width', 
#                      main = "Peak sharpness", prob = T) 
# lines(density( df_kvals$width), lwd = 2, col = "blue")

# Example -- Jan 2023
# 2. Dispersion stats 
# dispersion_stats <- parseDispersion(pathway_list_dispersion )
# 
# 3. Master data.frame 
# res = plot_dispersion_distribution(dispersion_stats , 
#                                    labels_idx = c(9,10,12,13,14,1,3), min_silh = 0.3, 
#                                    min_n_cluster = 2,  df_kvals = df_kvals , use_mean = F, 
#                                    include_metabolic = F, use_palette = 'Paired')

# 4. visualization
# res[[1]] + theme(text = element_text(size = 20)) + 
#     ylab('N clusters / n genes') + xlab('Dispersion')
# master_df <- res[[2]] # will now include all the summary statistics for all pathways 
# master_df <- master_df %>% mutate(short_label = str_trunc(pathway_name, 10) ) 
# We can no visualize: 

## Visualization of summary metrics -- Jan 2023
library(ggrepel)
fig5_scatter <- function(master_df = data.frame(), 
                        filter_width = T, 
                        filter_range = F, 
                        add_marginal = F , max_width = 0.4, 
                        max_range = 50 # absolute width of the peak 
                        ){

  

  if(filter_width)
    p <- master_df %>% dplyr::filter( width < max_width ) %>% 
        ggplot(aes(x = k_gene_ratio, y = median_disp, size = 1/width, label = short_label)) + 
        theme(text= element_text(size = 20 ) ) + theme_classic() + 
        geom_point() + theme(legend.position = 'left') + xlab('N cluster p. gene') + ylab('Dispersion')

  if(filter_range)
      p <- master_df %>% dplyr::filter( width_range < max_range ) %>% 
        ggplot(aes(x = k_gene_ratio, y = median_disp, size = 1/width, label = short_label)) + 
        theme(text= element_text(size = 20 ) ) + theme_classic() + 
        geom_point() + theme(legend.position = 'left') + xlab('N cluster p. gene') + ylab('Dispersion')

  # add labels 
  p <- p + geom_text_repel(size =3)

  if(add_marginal)
    ggMarginal(p, type = 'density', groupColour = F, groupFill = F)

  return(p)
}