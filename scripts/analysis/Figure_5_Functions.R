library(ggrepel)

## Functions:
perc_k_finder <- function(z_score = c(),
                          percentile = 0.9
){
  idx_max_z = which.max(z_score)
  max_z = max(z_score, na.rm=T)
  
  vals = rev(z_score[idx_max_z:length(z_score)])
  
  ret_idx = which.min(abs(vals-(max_z * percentile)))
  
  return(length(vals) - ret_idx + idx_max_z)
  
}

# Returns a data.frame of pathways and optimal number of clusters
process_silh <- function(save_file = F,
                         silh_files = c(),
                         silh_res_dir = "",
                         min_expr = 0.2){
  
  pathway_list = lapply(silh_files, 
                        FUN = function(x) strsplit(x, "_silh_plt.RDS")[[1]]) %>% unlist()
  
  # prepare the data frame for pathways with output 
  df_kvals = data.frame(name = pathway_list, 
                        k = rep(0, length(pathway_list)))
  
  # Read output from all pathways and run the pipeline using 
  # the optimal number of clusters found above 
  for (i in 1:length(silh_files)){
    
    silh_result = readRDS(paste(silh_res_dir, 
                                silh_files[[i]], 
                                sep=""))
    
    res = silhouette_zscore(silh_result = silh_result,
                            min_expr = min_expr,
                            x_offset = 6,
                            min.y = -0.05, 
                            max.y = 0.55,
                            k_max = silh_result[[2]]$k %>% max()
    )
    
    df_kvals$k[[i]] = perc_k_finder(z_score = res[[2]], percentile = 0.9)                    
  }
  
  if(save_file)
    
    saveRDS(df_kvals, paste(silh_res_dir, "df_kvals.RDS", sep=""))
  
  return(df_kvals)
  
}

# Peak score function 
plot_peakedness_score <- function(z_score_raw = c() , 
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
  plot(z_score_raw, 
       type = "l", 
       main = pathway_name, 
       xlab = "Number of clusters (k)")
  
  lines(z_score , 
        type = "l", 
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
  plot(abs(z_score - (max_z * percentile))/max_z, 
       type = "l", 
       xlab = "Number of clusters (k)")
  abline(h = peak_width_threshold, 
         col ="red")
  
  # absolute distance of z-score values relative to the peak 
  peak_dist <- abs(z_score - (max_z * percentile))/max_z
  
  # After applying the smoothing filter, we get some NA values 
  peak_dist_vals <- peak_dist[!is.na(peak_dist)]
  
  # how many k values have a z-score within 10% of the peak 
  # if too many k values have a z-score close to the max (peak) then we can classify 
  # this peak as "broad" 
  peak_width <- sum(peak_dist_vals <= peak_width_threshold)
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

peak_width_scores <- function(use_percentile = 0.95, # threshold for considering values within the peak range (% of max z-score)
                              smooth_window = 3,     # filtering window size for smoothing the z-score data before finding the max value. 
                              pct_opt = 0.9,         # will be deprecated!!: : percentile of max z-score that will be used to identify the optimal k - this can be different from the above percentile 
                              peak_threshold = 0.1,  # relative distance from the max to consider a k value within the peak 
                              silh_files = c(),
                              silh_res_dir = "",
                              min_expr = 0.2
){
  
  pathway_list = lapply(silh_files, 
                        FUN = function(x) strsplit(x, "_silh_plt.RDS")[[1]]) %>% unlist()
  
  # prepare the data frame for pathways with output 
  df_kvals = data.frame(name = pathway_list, 
                        k = rep(0, length(pathway_list)), 
                        width = rep(0, length(pathway_list)))
  
  # Read output from all pathwways and run the pipeline using 
  # the optimal number of clusters found above 
  for (i in 1:length(silh_files)){
    
    silh_result = readRDS(paste(silh_res_dir, 
                                silh_files[[i]], 
                                sep=""))
    
    res = silhouette_zscore(silh_result = silh_result,
                            min_expr = min_expr,
                            x_offset = 6,
                            min.y = -0.1, 
                            max.y = 0.6,
                            k_max = silh_result[[2]]$k %>% max()
    )
    
    pdf(paste(silh_res_dir, 
              pathway_list[[i]], 
              '_silhPeak.pdf', 
              sep = ""))
    
    # Computes the peak width 
    peak_width_list  <- plot_peakedness_score(z_score_raw = res[[2]], 
                                              pathway_name = pathway_list[[i]], 
                                              fil_win = smooth_window, 
                                              percentile = use_percentile, 
                                              peak_width_threshold = peak_threshold)
    
    dev.off() 
    df_kvals$k[[i]] = perc_k_finder(z_score = res[[2]], 
                                    percentile = pct_opt)  # 0.9 is the default value --- used in the manuscript 
    df_kvals$width[[i]] <- peak_width_list[[1]]
    df_kvals$min_k_peak[i] <- peak_width_list[[2]]
    df_kvals$max_k_peak[i] <- peak_width_list[[3]]
    
  }
  
  return(df_kvals) 
  
}

# 1. For each pathway, run the pipeline using the optimal number of clusters identified by process_silh() 
# 2. Export a data.frame of pathway, k, dispersion 
# 3. Save the full output as .RDS for each pathway 
computeDispersions <- function(use_min_k = F, # whether to choose the optimal k by considering the minimum value of k that is close to the peak (or the max -- default)
                               overwrite_files = F, #  if we are running only new pathways and we don't want to re-run pathways that finished already 
                               dispersion_dir = "",
                               silh_res_dir = "",
                               df_kvals = data.frame(),
                               silh_files = c()
){
  
  pathway_list = lapply(silh_files, 
                        FUN = function(x) strsplit(x, "_silh_plt.RDS")[[1]]) %>% unlist()
  
  # subset only pathways for which we don't have DISPERSION results already 
  if(!overwrite_files){
    
    #check which pathways have dispersion results already (if we don't want to overwrite data)
    done_pathways <- list.files(dispersion_dir)
    
    # exclude them from the list of pathways to run 
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
    
    saveRDS(control_res, 
            paste(dispersion_dir, 
                  pathway_list_dispersion[[i]], 
                  "_diversity.RDS", 
                  sep=""))
    
    print(paste("DONE ", 
                pathway_list_dispersion[[i]], 
                Sys.time()))
    
  }
  
}

# New code to compute dispersion 
parseDispersion <- function(pathway_list_dispersion = list(),
                            dispersion_dir = ""
){
  dispersion_list = list() 
  
  for(p in pathway_list_dispersion){
    
    path_disp = readRDS(paste0(dispersion_dir, 
                               p, 
                               '_diversity.RDS'))
    
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
    dplyr::filter(mean_silh > min_silh, n > min_n_cluster) %>% 
    group_by(pathway_name) %>% 
    summarise(silh = mean(mean_silh), 
              median_disp = median(diversity), 
              n_clusters = n())
  
  select_pathways = pathway_list$pathway[labels_idx]
  
  df_plot_labeled <- df_plot %>% 
    left_join(df_kvals %>% 
                rename(pathway_name = name), 
              by ="pathway_name") %>%
    mutate(pathway_label = ifelse(pathway_name %in% 
                                    select_pathways,
                                  str_trunc(pathway_name, 12) ,"")) %>% 
    dplyr::filter(k > 3)
  
  df_plot_labeled <- df_plot_labeled %>% 
    mutate(pathway_class = ifelse(grepl("bolis|ycle|Gluconeogenesis|ynthesis|Ammonia|Succinate|Fructose|Lysine|Pentose|Fatty Acids", 
                                        pathway_name),
                                  "metabolism","signaling/protein") )
  
  # LPA pathway appears 3 times 
  df_plot_labeled <- df_plot_labeled %>% dplyr::filter(!grepl(pattern = "LPA4|LPA5|LPA1|LPA3|LPA6", 
                                                              pathway_name)) 
  
  # merge with number of gnees 
  df_plot_labeled %>% 
    left_join(pathway_df %>% 
                rename(pathway_name = pathway) %>% 
                group_by(pathway_name) %>% 
                summarise(n_genes = n()), 
              by = "pathway_name") %>% 
    mutate(k_gene_ratio = k/n_genes) -> df_plot_labeled
  
  if(!include_metabolic)
    df_plot_labeled <- df_plot_labeled %>% 
    dplyr::filter(pathway_class != 'metabolism')
  
  g <- df_plot_labeled %>%
    ggplot(aes(x = median_disp, 
               y = k_gene_ratio, 
               size = silh, 
               label = pathway_label, 
               color = pathway_class)) + 
    geom_point() + 
    geom_text()
  
  return(list(g, df_plot_labeled))
}

plotly_dispersion_distribution <- function(dispersion_stats = data.frame(), 
                                           labels_idx = c(9,10, 12,13,14),
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
  
  select_pathways = pathway_list$pathway[labels_idx]
  
  df_plot_labeled <- df_plot %>% 
    left_join(df_kvals %>% 
                rename(pathway_name = name), 
              by ="pathway_name") %>%
    mutate(pathway_label = ifelse(pathway_name %in% select_pathways,
                                  str_trunc(pathway_name, 8) ,"")) %>% 
    dplyr::filter(k > 3)
  
  df_plot_labeled <- df_plot_labeled %>% 
    mutate(pathway_class = ifelse(grepl("bolis|ycle|Gluconeogenesis|Fructose|Lysine|Pentose|Fatty Acids", 
                                        pathway_name),
                                  "metabolism","signaling/protein"))
  
  # LPA pathway appears 3 times 
  df_plot_labeled <- df_plot_labeled %>% dplyr::filter(!grepl(pattern = "LPA2|LPA4", 
                                                              pathway_name)) 
  
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
                               line = list(width = 2)))
  return(list(fig, df_plot_labeled))
  
}

# Visualization of summary metrics -- Jan 2023
fig5_scatter <- function(master_df = data.frame(), 
                         filter_width = T, 
                         filter_range = F, 
                         add_marginal = F, 
                         max_width = 0.4, 
                         max_range = 50 # absolute width of the peak 
){
  if(filter_width)
    p <- master_df %>% 
      dplyr::filter(width < max_width) %>% 
      ggplot(aes(x = k_gene_ratio, 
                 y = median_disp, 
                 size = 1/width, 
                 label = short_label)) + 
      theme(text= element_text(size = 20)) +
      theme_classic() + 
      geom_point() + 
      theme(legend.position = 'left') + 
      xlab('N cluster p. gene') + 
      ylab('Dispersion')
  
  if(filter_range)
    p <- master_df %>% 
      dplyr::filter(width_range < max_rang) %>% 
      ggplot(aes(x = k_gene_ratio,
                 y = median_disp, 
                 size = 1/width, 
                 label = short_label)) + 
      theme(text= element_text(size = 20)) + 
      theme_classic() + 
      geom_point() + 
      theme(legend.position = 'left') + 
      xlab('Number of profile clusters per pathway gene') + 
      ylab('Dispersion')
  
  # add labels 
  p <- p + geom_text_repel(size =3)
  
  if(add_marginal)
    ggMarginal(p, type = 'density', groupColour = F, groupFill = F)
  
  return(p)
}