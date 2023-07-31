## Directories for the files
source("./scripts/analysis/imports_new.R")
output_dir <- "./data/processed_data/silhouette_res/silh_rds/"

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
pathway_df$pathway  <- pathway_df$pathway %>% str_replace('/', ' ')
pathway_df$pathway  <- pathway_df$pathway %>% str_replace('\\(', ' ')
pathway_df$pathway  <- pathway_df$pathway %>% str_replace('\\)', ' ')

# List of pathways to do silhouette analysis on
# Make list of pathways 
pathway_df %>% 
  dplyr::group_by(pathway) %>% 
  dplyr::count() %>% 
  dplyr::filter(n>7) %>% 
  dplyr::pull(pathway) -> filter_pathway_list 

# # Parameters to show
min_expr_threshold <- 0.3
min_genes_pathway <- 2

# Run pipeline in parallel
runPipeline <- function(pathway_name ="",
                        file_dir = ""){
  print(paste("Running ..", pathway_name , " " , Sys.time()))
    
  pathway_genes <- genesPathway(pathway_name = pathway_name,
                                pathway_df = pathway_df,
                                seurat_obj = master_seurat)
  
	# if it founds an error, it won't break but it will save the pathway as NA
	# we can then read the outputs and re-run the NA with a lower number of maxk 
	# the most likely cause of error is too many clusters 
	silh_plt = tryCatch({silhouettePlot(pathway_genes = pathway_genes, 
	                                    min_genes_on = min_genes_on, 
	                                    min_expr = min_expr, 
	                                    n_bootstraps = 5,
	                                    seurat_obj = master_seurat)
      }, 
      error = function(e) { print(e); return ("NA")
      }, 
      finally = { print("error handled") } )
    
  # save the file! This is what we want!! 
  if(silh_plt != "NA"){
    saveRDS(silh_plt, paste(file_dir, pathway_name, "_silh_plt_ale.RDS", sep = ""))
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


#######################
# Niv Update Nov 5 2022
# 0. Run the silhouette functions for all pathways and save results ... DONE 
# 1. Read RDS output
# 2. Find the optimal number of clusters  --- calls perc_k_finder() 
# 3. Run the pipeline and compute dispersion 
process_silh <- function(save_file = F,
                         silh_res_dir = ""){
  # returns a data.frame of pathways and optimal number of clusters 

  saved_files <- paste(silh_res_dir, 
                       list.files(path=silh_res_dir, pattern = "_ale.RDS"), 
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

  # prepare the data frame for pathways with output 
  df_kvals = data.frame(name = pathway_list_dispersion, 
                        k = rep(0, length(pathway_list_dispersion)))

  # Read output from all pathwways and run the pipeline using 
  # the optimal number of clusters found above 
  for (i in 1:length(pathway_list_dispersion)){
    silh_result = readRDS(paste(silh_res_dir, 
                                pathway_list_dispersion[[i]], 
                                "_silh_plt_ale.RDS", 
                                sep = ""))
    
    res = silhouette_zscore(silh_result = silh_result,
                            min_expr = min_expr_threshold,
                            x_offset = 6,
                            min.y = -0.05, 
                            max.y = 0.55,
                            k_max = silh_result[[2]]$k %>% max()
    )
    
    df_kvals$k[[i]] = perc_k_finder(z_score = res[[2]])                    
  }

  if(save_file)
    saveRDS(df_kvals, paste(silh_res_dir, "df_kvals.RDS", sep=""))

  return(df_kvals)
  
}

# Run the pipeline <- let's do this part in parallel 
### 

# 1. For each pathway, run the pipeline using the optimal number of clusters identified by process_silh() 
# 2. Export a data.frame of pathway, k , dispersion 
# 3. Save the full output as .RDS for each pathway 
computeDispersions <- function(use_min_k = F, # whether to choose the optimal k by considering the minimum value of k that is close to the peak (or the max -- default)
                               overwrite_files = F, #  if we are running only new pathways and we don't want to re-run pathways that finished already 
                               dispersion_dir = "",
                               silh_res_dir = "",
                               df_kvals = data.frame()
){

  # Find all the pathway results 
  saved_files <- paste(silh_res_dir, 
                       list.files(path=silh_res_dir, pattern = "_ale.RDS"),
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

  # load only the pathways that actually finished (silhoutte)
  pathway_list_dispersion <- full_pathway_list[saved_idx]

  # subset only pathways for which we don't have DISPERSION results already 
  if(!overwrite_files){

    #check which pathways have dispersion results already (if we don't want to overwrite data)
    done_pathways <- list.files(dispersion_dir)
    # exclude them from the list of pahtways to run 
    done_idx <- sapply(pathway_list_dispersion, 
                       function(x){sum(grepl(x, done_pathways))}) %>% 
      unlist %>% 
      as.logical
    
    pathway_list_dispersion <- pathway_list_dispersion[!done_idx]
  }

  pathway_list_dispersion <- pathway_list_dispersion[2:length(pathway_list_dispersion)]
  
  for (i in 1:length(pathway_list_dispersion)){

    if(!use_min_k){
      optimal_k_pathway = df_kvals %>% 
        dplyr::filter(name == pathway_list_dispersion[i]) %>% 
        dplyr::pull(k)
    }else{
      optimal_k_pathway = df_kvals %>% 
        dplyr::filter(name == pathway_list_dispersion[i]) %>% 
        dplyr::pull(min_k_peak)
    }

    print(paste("Running", 
                pathway_list_dispersion[i], 
                "k =", 
                optimal_k_pathway, 
                Sys.time(), sep = " "))
    
    pathway_genes = genesPathway(pathway_name = pathway_list_dispersion[i],
                                 pathway_df = pathway_df, 
                                 seurat_obj = master_seurat)
    
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
                                      dist_metric = "euclidean")
    
    saveRDS(control_res, paste(dispersion_dir, pathway_list_dispersion[[i]], "_diversity.RDS", sep=""))

    print(paste("DONE ", pathway_list_dispersion[[i]], Sys.time()))
  
    }

}

# New code to compute dispersion 
parseDispersion <- function(pathway_list_dispersion = list(),
                            dispersion_dir = ""
                            ){
  dispersion_list = list() 

  for(p in pathway_list_dispersion){
      
      path_disp = readRDS(paste0(dispersion_dir, p, '_diversity.RDS'))
      
      df_path = path_disp[[4]]$pathway
      df_path$pathway_name  = p
      dispersion_list[[p]] = df_path
      
      }

  dispersion_stats = do.call(rbind, dispersion_list)

  return(dispersion_stats)
  
}
# make statistics 
plot_dispersion_distribution <- function(dispersion_stats = data.frame(), 
                                         labels_idx = c(9,10,12,13,14),
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
    summarise(silh = mean(mean_silh), 
              median_disp = median(diversity), 
              n_clusters = n())

  select_pathways = pathway_list[labels_idx ]

  df_plot_labeled <- df_plot %>% 
    left_join(df_kvals %>% 
                rename(pathway_name = name), 
              by ="pathway_name") %>%
    mutate(pathway_label = ifelse(pathway_name %in% 
                                    select_pathways,str_trunc(pathway_name, 12) ,"")) %>% 
    dplyr::filter(k > 3)

  df_plot_labeled <- df_plot_labeled %>% 
                        mutate(pathway_class = ifelse(grepl("bolis|ycle|Gluconeogenesis|ynthesis|Ammonia|Succinate|Fructose|Lysine|Pentose|Fatty Acids", pathway_name),
                                         "metabolism","signaling/protein") )

  # LPA pathway appears 3 times 
  df_plot_labeled <- df_plot_labeled %>% dplyr::filter(!grepl(pattern = "LPA4|LPA5|LPA1|LPA3|LPA6", pathway_name)) 

  # merge with number of gnees 
  df_plot_labeled %>% 
    left_join(pathway_df %>% 
                rename(pathway_name = pathway) %>% 
                group_by(pathway_name) %>% 
                summarise(n_genes = n()) , 
              by ="pathway_name") %>% 
    mutate(k_gene_ratio = k/n_genes) -> df_plot_labeled

  if(!include_metabolic)
    df_plot_labeled <- df_plot_labeled %>% dplyr::filter(pathway_class != 'metabolism')

 
  g <- df_plot_labeled %>%
      ggplot(aes(x = median_disp, 
                 y = k_gene_ratio, 
                 size = silh, label = pathway_label, 
                 color = pathway_class)) + 
      geom_point() + 
      geom_text()

  return(list(g, df_plot_labeled))
}


plotly_dispersion_distribution <- function(dispersion_stats = data.frame(), 
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
    summarise(silh = mean(mean_silh), 
              median_disp = median(diversity), 
              n_clusters = n())


  select_pathways = pathway_list[labels_idx]


  df_plot_labeled <- df_plot %>% 
    left_join(df_kvals %>% 
                rename(pathway_name = name), 
              by ="pathway_name") %>%
    mutate(pathway_label = ifelse(pathway_name %in% select_pathways,
                                  str_trunc(pathway_name, 8) ,"")) %>% 
    dplyr::filter(k > 3)

  df_plot_labeled <- df_plot_labeled %>% 
                        mutate(pathway_class = ifelse(grepl("bolis|ycle|Gluconeogenesis|Fructose|Lysine|Pentose|Fatty Acids", pathway_name),
                                         "metabolism","signaling/protein"))

  # LPA pathway appears 3 times 
  df_plot_labeled <- df_plot_labeled %>% dplyr::filter(!grepl(pattern = "LPA2|LPA4", pathway_name)) 

  # merge with number of gnees 
  df_plot_labeled %>% 
    left_join(pathway_df %>% 
                rename(pathway_name = pathway) %>% 
                group_by(pathway_name) %>% 
                summarise(n_genes = n()), 
              by ="pathway_name") %>% 
    mutate(k_gene_ratio = k/n_genes) -> df_plot_labeled

  fig <- plot_ly(data = df_plot_labeled, 
                 x = ~median_disp, 
                 y = ~k_gene_ratio, 
                 color = ~pathway_class, 
                 colors = 'Set1',
                 text = df_plot_labeled$pathway_name,
                 hoverinfo = 'text',
                 marker = list(size = 10,
                               #color = 'rgba(255, 182, 193, .9)',
                               line = list(#color = 'rgba(152, 0, 0, .8)',
                                 width = 2)))
  return(list(fig, df_plot_labeled))

  }



#res = plot_dispersion_distribution(dispersion_stats , labels_idx = c(9,10,12,13,14,1,3), min_silh = 0.3, min_n_cluster = 2,  df_kvals = df_kvals , use_mean = F, include_metabolic = F, use_palette = 'Paired')
 
#res[[1]] + theme(text = element_text(size = 20)) + ylab('N clusters / n genes') + xlab('Dispersion')


# Compute peak width scores for all pathways 

# Peak score function 
plot_peakdness_score <- function(z_score_raw = c() , 
                                 pathway_name = "", 
                                 fil_win = 5, 
                                 percentile = 0.9,  # % of max to use as the peak height -- smaller values result in broader peaks 
                                 peak_width_threshold  = 0.1){
    
  par(mfrow = c(2,1))
    
  # smooth the z-score 
  z_score <- stats::filter(z_score_raw, 
                           rep(1,fil_win)/fil_win, 
                           sides = 1) 
  # absolute max 
  idx_max_z = which.max(z_score)
  max_z = max(z_score, na.rm=T)
  
  # Plot
  plot(z_score_raw, type = "l", 
       main = pathway_name)
    
  lines(z_score , type = "l", 
        col = "red")
  abline(v = idx_max_z)
  abline(h = max_z * percentile, 
         col ="blue")
  # plot the optimal k as detected by our percentile function 
  abline(v = perc_k_finder(z_score_raw, 
                           percentile),
         lty = 2, 
         col = "blue")  # only for visualization
  
  #subplot 2 -- relative difference to the peak percentile 
  plot(abs(z_score - (max_z * percentile))/max_z, type = "l")
  abline(h = peak_width_threshold  , col ="red")
  
  # absolute distance of z-score values relative to the peak 
  peak_dist <- abs(z_score - (max_z * percentile))/max_z
  
  # After applying the smoothing filter, we get some NA values 
  peak_dist_vals <- peak_dist[!is.na(peak_dist)]
  
  # how many k values have a z-score within 10% of the peak 
  # if too many k values have a z-score close to the max (peak) then we can classify 
  # this peak as "broad" 
  peak_width <- sum(peak_dist_vals <= peak_width_threshold )
  peak_width <- peak_width/length(peak_dist_vals) # this should account for pathways with k_max < 200 (default for most pathways)
  
  # let's save the range values 
  # min k value with high silhouette 
  min_k <- which(peak_dist_vals <= peak_width_threshold)[1] + fil_win
  # max value in the range 
  max_k <- rev(which(peak_dist_vals <= peak_width_threshold) )[1] + fil_win  
  
  # let's add vertical lines for the range of the peak 
  abline(v= min_k, col = 'blue')
  abline(v = max_k, col = 'blue') 
  
  return(list(peak_width, min_k, max_k)) # retuns a value after plotting 

  }

#peak_width_scores(silh_res_dir = "./data/processed_data/silhouette_res/silh_rds/")

# Will generate and save plots of peaks for all pathways 
# Returns a data.frame with the width scores and the optimal k for all pathways 
# NOTE: pathways with low peak scores should not be considered as their "optimal k" is not a good approximation 
peak_width_scores <- function(use_percentile = 0.95, # threshold for considering values within the peak range (% of max z-score)
                              smooth_window = 3,     # filtering window size for smoothing the z-score data before finding the max value. 
                              pct_opt = 0.9,         # will be depracted!!: : percentile of max z-score that will be used to identify the optimal k - this can be different from the above percentile 
                              peak_threshold = 0.1,   # relative distance from the max to consider a k value within the peak 
                              silh_res_dir = ""
                              ){
  saved_files <- paste(silh_res_dir, 
                       list.files(path=silh_res_dir, pattern = "_ale.RDS"), 
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

  # prepare the data frame for pathways with output 
  df_kvals = data.frame(name = pathway_list_dispersion, 
                        k = rep(0, length(pathway_list_dispersion)), 
                        width = rep(0, length(pathway_list_dispersion)))

    # Read output from all pathwways and run the pipeline using 
  # the optimal number of clusters found above 
  for (i in 1:length(pathway_list_dispersion)){
      silh_result = paste(silh_res_dir, 
                          pathway_list_dispersion[[i]], 
                          "_silh_plt_ale.RDS", sep = "")
      
      res = silhouette_zscore(silh_result = readRDS(silh_result),
                              min_expr = min_expr_threshold,
                              min.y = -0.05, 
                              max.y = 0.55
      )
      
      pdf(paste('./scripts/figures/peak_analysis/', 
                pathway_list_dispersion[[i]], 
                '_silhPeak.pdf', 
                sep = ""))
      
      # Computes the peak width 
      peak_width_list  <- plot_peakdness_score(z_score_raw = res[[2]], 
                                               pathway_name = pathway_list_dispersion[[i]], 
                                               fil_win = smooth_window, 
                                               percentile = use_percentile, 
                                               peak_width_threshold = peak_threshold)
      
      dev.off() 
      df_kvals$k[[i]] = perc_k_finder(z_score = res[[1]], 
                                      percentile = pct_opt)  # 0.9 is the default value --- used in the manuscript 
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