library(dplyr)
library(ggpubr)
library(ggplot2)
library(cluster)
library(Seurat)
library(pheatmap)
library(stylo)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(grid)
library(parallel)
library(gridExtra)
library(viridis)
library(superheat)
library(DIZutils)
library(stringr)
library(RColorBrewer)

## Creates a black color palette with n distinct colors
black_pal <- function(n){
  colfunc <- colorRampPalette(c("white", "black"))
  return(colfunc(n))
}

# Function to make categorical palette using glasbey_colors
makeQualitativePal <- function(
  n,               # Number of unique colors in palette
  rand_order = F,  # Whether or not we want to randomize colors assigned to labels on every run
  skip = 0,        # Whether to skip any colors (skip = 1 skips white)
  tail_colors = F,
  glasbey_use = T  # Use the Glasbey palette
){
  
  if(glasbey_use){
    col_vector=glasbey
    names(col_vector)<- NULL
    
  }else{
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) # only 70 qualitative colors
    
  }
  
  if(n <= length(col_vector)){
    if(rand_order == T){
      return(sample(col_vector, n))
      
    }else{
      # To add diversity we can get the last n colors of the array. Useful when plotting two pathways
      if(tail_colors){
        x_col = tail(col_vector, n)
        
      }else{
        x_col = col_vector[(1 + skip):(n + skip)]
        
      }
      return(x_col)
    }
    
  }else{
    return(c(col_vector,col_vector[1:(n-length(col_vector))]))
  }
}

## Function to return the genes in a pathway stores in "all_pathways" object
genesPathway <- function(which_pathway = 'Notch'){
  
  this_pathway = all_pathways %>% dplyr::filter(pathway == which_pathway) %>% pull(gene)
  this_pathway = this_pathway[which(this_pathway %in% row.names(master_seurat))]
  
  return(this_pathway)
}

## Permute the rows for each column of data.frame
randomize_columns <- function(df,           # data.frame of gene expression values (rows are cells, columns are genes)
                              this_pathway  # Pathway genes list
){
  for(i in 1:length(this_pathway))
    df[, this_pathway[i]] = sample(df[, this_pathway[i]])
  return(df)
  
}

## Retrieve the log-norm gene expression values and the annotations
makeMainDataFrame <- function(this_pathway, master_seurat = c()){
  
  pathway_matrix <- FetchData(master_seurat, this_pathway) # Fetch log norm counts from Seurat object
  devel_adult <- cbind(pathway_matrix, master_seurat@meta.data %>%
                         dplyr::select(Tissue, age, dataset, cell_ontology_class, Cell_class, cell_id)) # Append metadata annotations
  
  row.names(devel_adult)<-1:dim(devel_adult)[1] # Set the rownames
  devel_adult$global_cluster = 1:dim(devel_adult)[1] # here we haven't filtered anything. Comes directly from Seurat obj
  
  return(devel_adult)
}

## Returns MinMax normalized counts along with metadata at specified saturating quantile.
normalizedDevel <- function(
  master_seurat = c(),    # Which Seurat object to use
  this_pathway,           # List of pathway names
  sat_val = 0.99,         # Saturating quantile
  fill_zero_rows = F ,    # If a gene has all 0's, fill with very small number
  which_datasets = 'both' # Include developmental, adult, or both datasets
){
  
  # Get the log norm gene counts and annotations
  devel_adult <- makeMainDataFrame(this_pathway, master_seurat)
  
  if(!'cell_id' %in% names(master_seurat@meta.data))
    devel_adult %>% mutate(cell_id = paste(global_cluster, dataset, sep = "_")) -> devel_adult
  
  # Filter profiles based on specification of which_datasets
  if(which_datasets == 'devel'){
    devel_adult <- devel_adult %>% dplyr::filter(dataset %in% c("E6.5_8.5_Chan" ,"E6.5_8.5_Marioni","E5.5_Nowotschin" ,"E9.5_11.5_Tang" ,"Forelimb_E10.5_15.0") )
  }else if(which_datasets=='adult'){
    devel_adult <- devel_adult %>% dplyr::filter(!dataset %in% c("E6.5_8.5_Chan" ,"E6.5_8.5_Marioni","E5.5_Nowotschin" ,"E9.5_11.5_Tang" ,"Forelimb_E10.5_15.0") )
  }
  
  # Get dataframe of only gene expression values in the subset of cell states we want
  x = devel_adult[, this_pathway]
  
  # Compute the saturating values at a specified sat_val quantile
  max_sat_gene = apply(x, 2, quantile, sat_val) # starts from x
  
  # If there are genes that saturate to a value of 0, then replace the saturating value with the max. value of its expression
  if(sum(max_sat_gene == 0) > 0){
    max_val_gene = apply(x, 2, max) # starts from x
    max_sat_gene[max_sat_gene == 0] <- max_val_gene[max_sat_gene == 0]
  }
  
  # Actually saturate the values now
  for(s in 1:dim(x)[2])
    x[which(x[, s] > max_sat_gene[s]), s]<- max_sat_gene[s]
  
  # Apply MinMax Scaling to scale values between 0 to 1
  x <- min.maxNorm(x)
  
  # We can't do cosine distance on all zero genes so we just fill these with a very small value if we want
  if(fill_zero_rows)
    x[x == 0] = 10^-10
  
  # Modify the dataframe with the metadata to now have the MinMax scaled, saturated counts
  devel_adult[, this_pathway] <- x
  
  row.names(devel_adult) <- devel_adult$global_cluster # Reset column names
  
  return(devel_adult)
}

## MinMax Scaling of values
min.maxNorm <- function(x){
  maxs = apply(x, 2, max) # applies max operation to columns
  mins = apply(x, 2, min)
  
  # Scale the values
  for(i in 1:dim(x)[2]){
    x[, i] = (x[, i] - mins[i]) / (maxs[i] - mins[i])
  }
  
  return(x)
}

## quickPipeline processes the counts from an object with a specified optimal number of clusters
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
    this_pathway = all_pathways %>% dplyr::filter(pathway == which_pathway) %>% pull(gene)
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

## Build a heatmap with pathway expression profiles and class labels
quickHeatmap <- function(master_seurat = c(),
                         which_pathway = 'Bmp',       # Pathway name to draw from all_pathways object
                         filter_profiles = 'adult',   # Seurat object
                         k = 36,                      # k value for clusters
                         silent = F,                  #
                         min_genes_ = 1,              # Min. number of genes for pathway to be "ON"
                         min_expr_gene = 0.2,         # Min. expression cutoff for gene to be "ON"
                         save_pdf = F,                # save PDF
                         pdf_file = '',               # PDF file name (without directory)
                         return_heatmap = F,          # whether or not to plot
                         cluster_genes = T,           # Cluster genes
                         skip_plot = F                # Skip plotting but return the data.frame of values
){
  
  # Grab the pathway genes
  this_pathway = genesPathway(which_pathway)
  
  # Running the quickPipeline function to get the counts and the annotations
  res_list = quickPipeline(which_pathway = which_pathway,
                           master_seurat = master_seurat,
                           which_profiles = filter_profiles,
                           k_final = k,
                           manual_filter_cells = c(),
                           verbose = T,
                           min_genes_on = min_genes_,
                           min_expr = min_expr_gene)
  
  df_devel = res_list$df_devel
  df_devel$class_label <- as.character(df_devel$class_label)
  
  annoCols <- makeQualitativePal(k, glasbey_use = T, skip = 1)
  names(annoCols) <- unique(df_devel$class_label)
  annoCols <- list("class_label" = annoCols)
  
  # skip_plot will not plot the heatmap but just returns the data.frame
  if(!skip_plot){
    
    if(!save_pdf){
      dist_mat = dist.cosine(df_devel[, this_pathway] %>% as.matrix)
      dist_mat = dist(df_devel[, this_pathway])
      
      p = pheatmap(df_devel[, this_pathway],
                   annotation_row = df_devel %>% select(class_label),
                   annotation_colors = annoCols,
                   show_rownames = F,
                   fontsize = 10,
                   cutree_rows = k,
                   clustering_method = 'ward.D2',
                   clustering_distance_rows = dist_mat,
                   silent = silent,
                   col = black_pal(100),
                   cluster_cols = cluster_genes)
      
    }else{
      p = pheatmap(df_devel[,this_pathway ],
                   annotation_row = df_devel %>% select(class_label),
                   annotation_colors = annoCols,
                   show_rownames = F, fontsize = 10,
                   cutree_rows = k, clustering_method = 'ward.D2',
                   clustering_distance_rows = dist.cosine(df_devel[,this_pathway ] %>% as.matrix ),
                   silent = silent,
                   col = black_pal(100),
                   cluster_cols = cluster_genes,
                   filename = pdf_file,
                   height = 5,
                   width = 5)
      
    }
    if(return_heatmap){
      return(p)
    }else{
      return(df_devel)
    }
  }else{
    return(df_devel)
  }
}

## Takes an array of different bootstrap runs and returns the average silhouette score for each k value from all runs
makeBootstrap_df <- function(s_boots){
  
  # rbind, so each column represents k value
  boot_mat  = do.call(rbind, s_boots)
  
  boot_df = data.frame(m = apply(boot_mat, 2, mean), s = apply(boot_mat, 2, sd)) # Mean over columns
  boot_df$k = 1:dim(boot_df)[1]
  
  return(boot_df)
}

## Bootstraps silhouette scores from from global expression profiles
silhPathwayBootstrap <- function(
  x,                         # Pathway name
  master_obj = c(),          # Seurat object
  cells_include = c(),
  k_max = 100,               # Maximum values of k_cutoffs to compute the silhouette score for
  clust_method = "complete", # Clustering clust_method
  sat_val = 0.99,            # Saturating quantile
  dist_metric = 'euclidean', # Distance metric
  pct_boots = 0.9,           # Fraction of cells we want to sample from dataset
  n_boots = 100,             # Number of bootstraps to run
  weighted_silhouette = F,   # Whether to weight silhouette score by the number of rows in each cluster
  control_silh = F           # Whether this is a bootstrapping for clusters or a negative control that randomzes the data
){
  this_pathway = x # Get pathway genes
  devel_adult_main <- normalizedDevel(this_pathway,
                                      sat_val,
                                      fill_zero_rows = T,
                                      master_seurat = master_obj)
  
  boot_list = list() # Bootstrap values list
  for(b in 1:n_boots){
    
    if(!control_silh){
      # Running bootstraps on the pathway genes
      # Filter out cell types that are not expressing.
      # but sample the list of cell types before to create a bootstrap distribution
      cells_include_boots  = sample(cells_include, round(length(cells_include) * pct_boots))
      
      devel_adult_main %>% dplyr::filter(cell_id %in% cells_include_boots) -> devel_adult
      row.names(devel_adult) <- devel_adult$global_cluster
      
    }else{
      # Control bootstrap: re-shuffle the data to get a lower bound on silhouette score
      # by destroying gene-gene correlations but keeping the distribution of each gene.
      # Before, we still need to filter for the cell types to include
      devel_adult_main %>% dplyr::filter(cell_id %in% cells_include) -> devel_adult_main
      
      # Randomize expression within those cells! We don't want gene expression values from cell types not originally included
      devel_adult = randomize_columns(devel_adult_main, this_pathway)
      
    }
    
    x = devel_adult[, this_pathway] %>% as.matrix
    
    # Compute distance metric on gene expression values
    if(dist_metric =="euclidean"){
      dist_mat = dist(x)
      
    }else if(dist_metric=="cosine"){
      dist_mat = dist.cosine(x) # Matrix should not have zero rows
      
    }
    
    # We compute the clustering only once
    main_clustering = hclust(dist_mat, method = clust_method)
    ss = c() # Stores the silhouette scores
    for(k in 2:k_max){
      # Cut the tree for each k value
      clustering = cutree(main_clustering, k)
      
      devel_adult$class_label = clustering
      
      master_clustered = devel_adult # data.frame with data and motif labels
      
      s1 <- cluster_silhouette(master_clustered,
                               this_pathway = this_pathway, # Pathway genes
                               dist = dist_metric # original clustering
      )
      
      # Two ways to compute the silhouette score to return
      if(!weighted_silhouette){
        # normal silhouette score -- average across all clusters
        ss[k] = mean(s1$ms)
        
      }else{
        # A weighted average to consider the number of cell types within each cluster
        ss[k] = mean(s1$ms * s1$n)
        
      }
    }
    
    boot_list[[b]] = ss
  }
  
  return(boot_list)
}

## Function that computes silhouette scores for a given number of clusters
cluster_silhouette <- function(df = data.frame(), # data.frame in wide format with genes and class labels
                               this_pathway = c(),# pathway genes
                               dist = 'euclidean',
                               return_singles = F # directly return the mean silhouette score
){
  
  # Assumes that the data.frame contains the pathway profiles in wide format
  # The input data.frame must contain global_cluster AND class_label
  ## class_label must be numeric friendly
  df %>% dplyr::select(c(this_pathway, global_cluster, class_label)) -> df_mat
  
  # Silhouette requires numeric cluster labels
  labs <- df_mat$class_label %>% as.numeric()
  names(labs) <- df_mat$global_cluster
  x <- df_mat[, this_pathway] %>% as.matrix
  row.names(x) <- df_mat$global_cluster
  
  # Computes the silhouette scores using specified distance metric
  if(dist == 'cosine'){
    s <- silhouette(labs, dist.cosine(x))
    
  }else if(dist == 'euclidean'){
    s <- silhouette(labs, dist(x))
    
  }
  
  s_df <- data.frame(class_label = s[, 1], silh = s[, 3])
  
  s_df_mean <- s_df %>% dplyr::group_by(class_label) %>%
    dplyr::summarise(ms = mean(silh)) %>%
    arrange(desc(ms)) %>% as.data.frame() %>%
    left_join(df %>% dplyr::group_by(class_label) %>% count, by="class_label")
  
  if(!return_singles){
    return(s_df_mean)
    
  }else{
    df$silh = s_df$silh
    return(df)
    
  }
}

## Compute pathway gene and randomized silhouette scores and plots them together
silhouettePlot <- function(which_pathway = 'Notch',   # Specify pathway name
                           min_ON = 2,                # # Min. number of genes for pathway to be "ON"
                           min_expr = 0.25,           # Min. expression cutoff for gene to be "ON"
                           n_bootstraps=100,          # Number of bootstrap replicates
                           pct_bootstrap = 0.9,       # Percent of cells to ??
                           max_k = 100,               # Max. number of clusters to compute silhouette scores for
                           clust_metric = "cosine",   # Clustering metric
                           clust_method = "ward.D2"   # Clustering method
){
  
  this_pathway = genesPathway(which_pathway); # Grab pathway genes
  
  # Get pathway MinMax normalized counts
  devel_adult = quickHeatmap(which_pathway = which_pathway,
                             master_seurat = master_seurat,
                             filter_profiles = 'both',
                             k = 20,
                             min_genes_ = min_ON,
                             min_expr_gene = min_expr)
  
  # 1. Compute the real pathway silhouette scores
  s_boots = silhPathwayBootstrap(this_pathway,
                                 devel_adult$cell_id,
                                 k_max= max_k,
                                 master_obj = master_seurat,
                                 dist_metric = clust_metric ,
                                 clust_method = clust_method,
                                 control_silh = F,
                                 n_boots = n_bootstraps,
                                 pct_boots = pct_bootstrap)
  
  boot_df = makeBootstrap_df(s_boots)
  
  # 2. Compute the negative control silhouette scores
  s_boots = silhPathwayBootstrap(this_pathway,
                                 devel_adult$cell_id,
                                 k_max= max_k,
                                 master_obj = master_seurat,
                                 dist_metric = clust_metric,
                                 clust_method = clust_method,
                                 control_silh = T,
                                 n_boots = n_bootstraps)
  
  boot_df_control = makeBootstrap_df(s_boots)
  
  # 3. Plot together
  boot_df %>% ggplot(aes(x = k, y = m)) +
    geom_ribbon(aes(ymin = m - s, ymax = m + s), fill = 'lightblue', alpha = 0.2) +
    geom_line(color ='blue') + geom_ribbon(data = boot_df_control, aes(ymin = m - s, ymax = m + s), alpha = 0.2) +
    geom_line(data = boot_df_control)  + theme_pubr(base_size = 12) + ylab('Silhouette') +
    xlab('Number clusters') + ggtitle(which_pathway) -> g
  
  # Returns plot and both data frames
  return(list(g, boot_df, boot_df_control))
}

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
    geom_vline(xintercept = which.max(z_score), linetype = "dashed") + theme_pubr(base_size = 10) +
    ylab("Z-score") + xlab("Number clusters") + theme(axis.title.y = element_text(colour = , color = rgb(54, 55, 149, max = 255)))
  
  # Merge the two plots into 1 using ggarrange
  g <- ggarrange(p1, NULL, p2, ncol = 1, align = "hv", heights = c(1, -0.175, 0.5))
  
  return(g)
}

## Compute silhouette score for given labels, and returns mean silhouette score per class
clusteringSilhouette <- function(df = data.frame(),     # data.frame in wide format with genes and class labels
                                 dist_metric = 'cosine',# distance metric
                                 gene_list = c()        # List of genes in pathway
){
  # Compute the silhouette score for k cluster labels
  if(dist_metric == 'cosine'){
    d = silhouette(df$class_label %>% as.numeric, dist = dist.cosine(df[, gene_list] %>% as.matrix))
    
  }else{
    d = silhouette(df$class_label %>% as.numeric, dist = dist(df[, gene_list] %>% as.matrix))
    
  }
  
  # Store silhouette results in data.frame of class label, neigbhoring cluster, and silhouette score
  silh_df = data.frame(class_label = d[, 1], neighbor = d[, 2], silh = d[, 3])
  
  # Take the silhouette score for a cluster as the average of that of all its members
  silh_df %>% group_by(class_label) %>% summarise(mean_silh = mean(silh)) -> silh_df
  
  return(silh_df)
}

## Computes pathway and global diversity, along with mean ss for each class
globalClustering <- function(df_devel = data.frame(),   # data.frame in wide format with genes and class labels
                             master_seurat = c(),       # Seurat object
                             k_final = 25,              # k value for clusters
                             n_pcs = 100,               # Number of PCs to use
                             manual_embedding = c(),    # Matrix should represent a latent space (i.e. PCA, UMAP), do not provide if using n_pcs
                             dist_metric = 'euclidean', # Distance metric
                             pathway_genes = c()        # List of pathway genes
){
  
  df_devel$class_label <- df_devel$class_label %>% as.character()
  
  # Rename df_devel by cell_id
  row.names(df_devel) <- df_devel$cell_id
  
  # Set class_label colors
  colors_1206$class_label = makeQualitativePal(length(df_devel$class_label %>% unique))
  names(colors_1206$class_label) <- df_devel$class_label %>% unique
  
  # Get the Embedding
  # umap_coords by default uses cell_id as the row.name
  if(length(manual_embedding) == 0){
    # If user specifies n_pcs
    umap_coords <- Embeddings(master_seurat, reduction = 'pca')
    umap_coords <- umap_coords[, 1:n_pcs]
    
  }else{
    # Else use the user's Embedding, when provided directly, n_pcs is ignored
    # row.names equal cell_id
    umap_coords = manual_embedding
    
  }
  
  # Use the cells that only are in our data.frame object, with the pathway "ON"
  scaled_data = umap_coords[df_devel$cell_id, ]
  
  # Compute distance in embedding space between cells
  if(dist_metric =='euclidean')  user_dist = dist(scaled_data)  else  user_dist = dist.cosine(scaled_data)
  
  # Global clustering of cell types uses euclidean distance + ward.D2
  p_global = pheatmap(scaled_data %>% t ,
                      annotation_col = df_devel %>% dplyr::select(age, dataset, Tissue, Cell_class, class_label),
                      annotation_colors = colors_1206,
                      clustering_method = 'ward.D2',
                      cutree_cols = k_final,
                      show_rownames = F,
                      show_colnames = F,
                      fontsize = 12,
                      silent = T,
                      clustering_distance_cols = user_dist)
  
  # For each pathway label, we calculate the mean distance in the embedding space for the cell types in that label
  
  k_final = length(df_devel$class_label %>% unique)
  
  # Distance metric to be used as matrix
  # We convert the distance to a matrix so we can compute the pair-wise similarity
  dist_mat = as.matrix(user_dist)
  
  # For each class, compute the pairwise similarity
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$class_label == x]
    
    if(length(class_cells) > 1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% max())
      
    }else{      # If there is only one member in the class, then pairwise similarity is zero
      return(0)
    }
    
  }) -> diversity_pathway
  
  # Perform the same pairwise similarity computation on global labels
  
  global_labels = cutree(p_global$tree_col, k = k_final)
  df_devel$global_label = global_labels
  
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$global_label == x]
    
    if(length(class_cells)>1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% max())
      
    }else{
      return(0)
    }
    
  }) -> diversity_global
  
  # Transcriptome diversity
  df_devel %>% group_by(global_label) %>% count %>% as.data.frame() -> global_stats
  global_stats$diversity = diversity_global
  
  # Pathway diversity
  df_devel %>% group_by(class_label) %>% count %>% as.data.frame() %>%
    mutate(class_label = as.numeric(class_label)) %>%
    arrange(class_label) -> path_stats
  path_stats$diversity = diversity_pathway
  
  # Compute silhouette scores for this clustering
  # By default, we perform gene clustering using the cosine distance and ward.D2
  # so we also use cosine distance for the silhouette score
  
  pathway_silh = clusteringSilhouette(df = df_devel, dist_metric = 'cosine', gene_list = pathway_genes)
  
  # Join with the main data.frame
  path_stats %>% left_join(pathway_silh, by = 'class_label') -> path_stats
  # Dummy column to match later
  global_stats$mean_silh = 0
  
  path_stats$type = "pathway"
  global_stats$type = "transcriptome"
  
  path_stats %>% rename(label = class_label) -> path_stats
  global_stats %>% rename(label = global_label) -> global_stats
  
  # Rank the profiles by diversity
  path_stats$rank = rank(path_stats$diversity)
  
  return(list('pathway' = path_stats, 'global' = global_stats))
}

## Positive control that scrambles profiles to compute what recurrent dispersion would look like
positive_control <- function(df_devel = data.frame(),
                             master_seurat = c(),
                             n_pcs = 100,
                             n_random = 100,
                             user_embedding = c(),
                             dist_metric  = 'euclidean'
){
  # Cells are filtered
  df_devel$class_label <- df_devel$class_label %>% as.character()
  k_final = length(df_devel$class_label %>% unique)
  
  # PCA embedding
  # Number of principal components
  # umap_coords by default uses cell_id as the row.name
  if(length(user_embedding) == 0){
    umap_coords <- Embeddings(master_seurat,reduction = 'pca')
    umap_coords <- umap_coords[, 1:n_pcs]
    
  }else{
    # Here the latent space already has the row.names
    # dimensions should be specified before
    umap_coords = user_embedding
  }
  
  # After we specified the embedding we are using for the global clustering
  # we set up the distance metric to be used: default euclidean
  if(dist_metric =='euclidean') user_dist = dist(umap_coords) else  user_dist = dist.cosine(umap_coords )
  
  # Convert to matrix. Keeps the same row.names in the original embedding
  # Just make sure the embedding is correctly named by rows
  dist_mat = as.matrix(user_dist)
  
  control_list = list()
  
  for(i in 1:n_random){
    df_devel$class_label = sample(df_devel$class_label) # Scramble the pathway labels
    
    # For each pathway label, we calculate the mean distance
    # in the PCA space for the cell types in that label
    sapply(1:k_final, function(x){
      
      class_cells = df_devel$cell_id[df_devel$class_label == x]
      
      if(length(class_cells) > 1){
        dist_mat[class_cells, class_cells] -> dist_class
        
        return(dist_class[upper.tri(dist_class)] %>% max())
        
      }else{
        return(0)
      }
      
    }) -> diversity_pathway
    
    # Pathway
    df_devel %>% group_by(class_label) %>% count %>% as.data.frame() %>%
      mutate(class_label = as.numeric(class_label)) %>%
      arrange(class_label) -> path_stats
    path_stats$diversity = diversity_pathway
    
    path_stats$type = "pathway"
    path_stats %>% rename(label = class_label) -> path_stats
    # rank the profiles by diversity
    path_stats$rank = rank(path_stats$diversity)
    
    control_list[[i]] = path_stats
  }
  return(control_list)
}

## Compute pathway and global ss
parallel_pipeline <- function(x,                          # data.frame of pathway expression values
                              cut_k = c(),                # k value for clusters
                              seurat_master = c(),        # Seurat object
                              profile_filter = 'adult',   # Include developmental, adult, or both datasets
                              manual_cell_types = c(),    # Any cell types to filter out
                              global_npcs = 100,          # Number of PCs to use if using PC embedding
                              user_embedding = c(),       # Else user-specified Embedding matrix
                              global_metric = 'euclidean' # Distance metric to use to compute similarities
){
  
  # quickPipeline takes manual_cell_types as the list with cell_ids to include
  # If empty, the function will filter cell_types based on the expression of the selected genes
  res_list = quickPipeline(which_pathway = x,
                           master_seurat = seurat_master,
                           which_profiles = profile_filter,
                           k_final = cut_k,
                           manual_filter_cells = manual_cell_types,
                           rand_control = T)
  
  stats_list = globalClustering(df_devel = res_list$df_devel,
                                master_seurat = seurat_master,
                                k_final = cut_k,
                                n_pcs = global_npcs,
                                manual_embedding = user_embedding,
                                dist_metric = global_metric,
                                pathway_genes = x)
  return(stats_list)
}

# Runs negative control tests for cell-type specific profiles by computing
# __ and diversity from non-pathway genes
controlDiversity <- function(null_dist = c(),
                             which_pathway = c(),
                             which_k = c(),
                             seurat_obj= c() ,
                             n_rand=100,
                             filter_profile = 'adult',
                             filter_manual_cells = c(),
                             verbose = F, n_pcs = 100,
                             manual_embedding = c(),
                             global_dist_metric = 'euclidean'
){
  
  control_stats = list()
  
  this_pathway = genesPathway(which_pathway)
  
  rand_list = list()
  for(rr in 1:n_rand)
    rand_list[[rr]] = sample(null_dist, length(this_pathway))
  
  # Run control pipeline in parallel for the whole list of random pathways
  # specify the manual list of cell types: overrides the filtering by min.expression
  control_stats = mclapply(rand_list, parallel_pipeline,
                           cut_k = which_k,
                           seurat_master = seurat_obj,
                           manual_cell_types = filter_manual_cells,
                           profile_filter = filter_profile,
                           global_npcs = n_pcs,
                           user_embedding = manual_embedding,
                           global_metric = global_dist_metric,
                           mc.cores = 10)
  
  # make control_df
  print(control_stats[[1]])
  do.call(rbind, lapply(control_stats, function(x){rbind(x$pathway %>% select(-rank), x$global)})) -> control_df
  
  return(control_df)
  
}


##
fullControlPathway <- function(this_pathway = c(),      # List of pathway genes
                               k_pathway = c(),         # k value for clusters
                               filter_pathway = 'adult',# Which profiles to keep
                               this_seurat = c(),       # Seurat object
                               null_list = c(),         # List of genes for control, HVGs to generate random data
                               n_samples = 100,
                               filter_manual = F,
                               verbose_output = F,
                               min_expr_gene = 0.2,     # Min. expression cutoff for gene to be "ON"
                               min_genes_ON = 1,        # Min. number of genes for pathway to be "ON"
                               n_pcs_global = 100,      # Number of PCs to use
                               embedding_matrix = c(),  # Number of principal components to consider for diversity metric
                               celltype_dist_metric = 'euclidean' # distance metric to be use in the global embedding (e.g., PCA) for global similarity across cell types
){
  
  # 1. Run the quick pipeline for the real pathway -- does not use any global metric for transcriptome
  #    with the specified k (found previously as optimal)
  res_list = quickPipeline(which_pathway = this_pathway,
                           master_seurat = this_seurat,
                           which_profiles = filter_pathway,
                           k_final = k_pathway, verbose = verbose_output,
                           min_genes_on = min_genes_ON, min_expr = min_expr_gene)
  if(verbose_output)
    print(paste(this_pathway,": ", dim(res_list$df_devel)[1]," expressing cells"))
  
  
  # 2. Compute diversity for the cell types containing each of the pathway profiles
  # Here we use an embedding to compute the global metric (PCA, VAE, ICA, etc)
  stats_list = globalClustering(df_devel = res_list$df_devel,
                                master_seurat = this_seurat,
                                k_final = k_pathway,
                                n_pcs = n_pcs_global,
                                manual_embedding = embedding_matrix,
                                dist_metric = celltype_dist_metric,
                                pathway_genes = genesPathway(this_pathway))
  
  
  # 3. Negative control: random sets of genes from null list
  # Note: here we want two types of filtering:
  #       current: filter cell types based on the expression of the random set
  #       probably better: filter cell types that express the real pathway so the diversity is comparable
  # new option to filter the same cell types as the real pathway
  if(filter_manual){
    manual_cell_types = res_list$df_devel$cell_id
  }else{
    manual_cell_types = c()
  }
  
  # Add the filter using filter_manual_cells parameter
  control_df = controlDiversity(null_dist = null_list,
                                which_pathway = this_pathway,
                                which_k = k_pathway,
                                seurat_obj = this_seurat,
                                n_rand = n_samples,
                                filter_profile = filter_pathway,
                                filter_manual_cells = manual_cell_types,
                                verbose = verbose_output,
                                n_pcs = n_pcs_global,
                                manual_embedding = embedding_matrix,
                                global_dist_metric = celltype_dist_metric)
  
  # Split negative control results: pathway (random genes), transcriptome (global diversity)
  control_df_path = control_df %>% dplyr::filter(type == 'pathway')
  control_df_global = control_df %>% dplyr::filter(type == 'transcriptome')
  
  
  # 4. Positive control: re-distribute pathway profiles in randomly selected cell types.
  # returns diversity scores only for the fake pathway (no global calculation)
  pos_control = positive_control(res_list$df_devel,
                                 master_seurat = this_seurat,
                                 n_random = n_samples,
                                 n_pcs = n_pcs_global,
                                 user_embedding = embedding_matrix,
                                 dist_metric = celltype_dist_metric
  )
  # Merge all in data.frame
  pos_control = do.call(rbind, pos_control)
  
  
  df_diversity = rbind(data.frame(d = control_df_path$diversity, type = 'null dist'), #random sets
                       data.frame(d = control_df_global$diversity, type = 'transcriptome dist'), #global dendrogram using same cells from random sets
                       data.frame(d = stats_list$pathway$diversity, type = 'pathway'), # actual pathway
                       data.frame(d = stats_list$global$diversity, type = 'transcriptome'), # global dendrogram using the same cells than the real pathway
                       data.frame(d = pos_control$diversity, type = 'pos control')) #randomized pathway profiles across cell types
  
  df_recurrence = rbind(data.frame(d = control_df_path$diversity * control_df_path$n, type = 'null dist'),
                        data.frame(d = control_df_global$diversity * control_df_global$n, type = 'transcriptome dist'),
                        data.frame(d = stats_list$pathway$diversity * stats_list$pathway$n, type = 'pathway'),
                        data.frame(d = stats_list$global$diversity * stats_list$global$n, type = 'transcriptome'),
                        data.frame(d = pos_control$diversity * pos_control$n, type ='pos control'))
  
  # We can also return the main data.frame, with class_label, diversity and rank for each profile
  clustered_data <- res_list$df_devel %>% left_join(stats_list$pathway %>%
                                                      rename(class_label = label) %>%
                                                      select(rank, diversity, class_label, n), by = 'class_label')
  
  return(list(diversity = df_diversity, recurrence = df_recurrence, profiles = clustered_data, rank = stats_list))
}

## Violin diversity plot
violin_diversity <- function(control_res = list(),
                             p_value = T,
                             plot_null = F
){
  
  my_comparisons = list(c('pathway','transcriptome'), c('pathway','pos control') )
  
  if(!plot_null){
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway', 'transcriptome', 'pos control'))
    lines_palette = c("#00AFBB", "#a8a8a8", "#000000") # color order: pathway, transcriptome, random
    
  }else{
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway', 'transcriptome', 'pos control', 'null dist'))
    lines_palette = c("#00AFBB", "#a8a8a8", "#000000",'#6f7e96') # color order: pathway, transcriptome, random
    
  }
  
  # Plot the diversity
  violin_plot <- ggviolin(df,
                          x = 'type',
                          y = 'd',
                          fill = 'type',
                          palette = lines_palette,
                          add = 'boxplot',
                          add.params = list(fill ='white'),
                          font.label = list(size = 10, size = "plain"))
  
  # Add p-value
  if(p_value){
    violin_plot + stat_compare_means(comparisons = my_comparisons) +
      stat_compare_means(label.y = control_res$diversity$d %>% max ) -> violin_plot
    
  }else{
    violin_plot + scale_y_continuous(trans='sqrt') + scale_x_discrete(labels = NULL) + coord_flip() -> violin_plot
    
  }
  
  violin_plot + ylab('Cell type diversity') +
    xlab('Cell type clustering') + theme_pubr(base_size = 10 ) -> violin_plot
  
  return(violin_plot)
}

## ECDF plotting function
ecdf_diversity <- function(control_res = list(),
                           plot_null = F,
                           pathway_name = 'Bmp'
){
  
  if(!plot_null){
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control'))
    
    lines_palette = c("#00AFBB", "#000000", "#a8a8a8") # color order: pathway, transcriptome, random
    
  }else{
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control','null dist'))
    lines_palette = c("#00AFBB", "#000000", "#a8a8a8",'#6f7e96') # color order: pathway, transcriptome, random
    
  }
  quantile_line = quantile(df %>% dplyr::filter(type=='transcriptome') %>% pull(d), 0.9)
  
  df  %>% ggplot(aes(x = d, col = type)) + stat_ecdf()  + ggtitle(pathway_name) +
    scale_x_continuous(trans = 'sqrt') + theme_pubr(base_size = 10) +
    scale_color_manual(values = lines_palette) +
    ylab('Fraction of clusters') + xlab('Cell state dispersion') +
    geom_hline(yintercept = 0.95,linetype="dashed", color = "lightgray", size = 1) -> ecdf_plot
}

## Plotting function for ranked diversity values
diversity_plot <- function(stats_list = list(), # Pathway and global diversity statistics from running fullControlPathway
                           title_main = 'Notch' # Pathway name for plot title
){
  # Extract pathway and global diversity values sorted by rank
  path_stats = stats_list$pathway
  global_stats = stats_list$global
  upper_quantile = 0.75
  
  g <- path_stats %>% ggplot(aes(x = rank, y = diversity)) + geom_point() +
    geom_hline(yintercept = quantile(global_stats$diversity, upper_quantile), linetype = "dashed", color = "lightgray", size = 1) +
    geom_hline(yintercept = quantile(global_stats$diversity, 1 - upper_quantile), linetype = "dashed", color = "lightgray", size = 1) +
    ylab("Cell type diversity (PCA)")  + theme_pubr(base_size = 10) +
    xlab("Pathway profile")  + ggtitle(title_main) +
    annotate(geom = "text", x = 3, y = quantile(global_stats$diversity, upper_quantile) + 2, label = "Expected", color="black")
  
  return(g)
}

## Actually compute the diversity values of the specified number of profiles and plot them
rank_diversity <- function(which_pathway ='Bmp',
                           k_motifs = 20,
                           min_expression = 0.2,
                           min_genes_pathway = 2,
                           embedding_matrix = c(),
                           global_dist_metric = 'euclidean',
                           make_plot = T
){
  
  # takes the pathway name as input:
  # all parameters come from upstream functions
  res_list = quickPipeline(which_pathway = which_pathway,
                           master_seurat = master_seurat,
                           k_final = k_motifs,
                           min_expr = min_expression,
                           min_genes_on = min_genes_pathway,
                           which_profiles = 'both')
  
  # Get the diversity values
  stats_list = globalClustering(df_devel = res_list$df_devel,
                                master_seurat = master_seurat,
                                k_final = k_motifs,
                                manual_embedding = embedding_matrix,
                                dist_metric = global_dist_metric,
                                pathway_genes = genesPathway(which_pathway))
  
  if(make_plot){
    return(diversity_plot(stats_list, which_pathway))
  }else{
    return(stats_list)
  }
}

## Computes the global transcriptome distance using a specified embedding and plots pathway class labels on global dendrogram
globalDendrogramPlot2 <- function(control_res = list(),                         # Main results object
                                  seurat_obj = c(),                             # Seurat object
                                  hvg_genes = c(),                              # List of genes to create the global dendrogram
                                  dist_metric ='cosine',                        # Clustering distance in global space. applies for expression and PCA spaces
                                  clust_method = "ward.d2",                     # Hierarchical clustering method
                                  use_pca = F,                                  # Whether to use PCA coordinates instead of gene expression for global dendrogram
                                  npcs = 1:100,                                   # If use_pca, then n_pcs to use
                                  save_pdf = T,                                 # Save plot
                                  save_dir = 'plots/figuresApril/fig4/'         # Directory to save plot to
){
  
  # Get the pathway states of the cells which have the pathway "ON"
  pathway_df = control_res$profiles
  row.names(pathway_df) <- pathway_df$cell_id
  
  meta_full <- master_seurat@meta.data %>% select(cell_ontology_class,
                                                  Cell_class, dataset, age)
  
  meta_full$cell_id <- row.names(meta_full)
  
  # Merge the pathway states with the whole dataset
  meta_full$class_label = 0 # assign cells with pathway "OFF" a class_label = 0
  meta_full[row.names(pathway_df), ]$class_label <- pathway_df$class_label
  
  pathway_df <- meta_full
  
  # Compute the global distance in the specified Embedding space
  if(!use_pca){
    # Use highly variable genes
    hvg_df = makeMainDataFrame(hvg_genes, master_seurat = seurat_obj)
    hvg_genes = hvg_genes[hvg_genes %in% colnames(hvg_df)]
    
    # Filter only for expressing cell types in this pathway
    hvg_df %>% dplyr::filter(cell_id %in% pathway_df$cell_id) -> hvg_df
    row.names(hvg_df) <- hvg_df$cell_id
    
    # Scale the gene expression values so mean is zero
    x = hvg_df[, hvg_genes]
    x_scaled = scale(x) # standardScaler
    row.names(x_scaled) <- hvg_df$cell_id
    
  }else{
    # Get the values in PCA space
    x_scaled = Embeddings(seurat_obj, reduction='pca')
    x_scaled = x_scaled[pathway_df$cell_id, npcs] # already set row names as cell id
    
  }
  
  # Compute distance in the global transcriptome space
  if(dist_metric == 'cosine') clust_dist_metric = dist.cosine(x_scaled) else clust_dist_metric = dist(x_scaled)
  
  k_motifs = length(pathway_df$class_label %>% unique) # number of pathway states
  
  # Colors for pathway classes --- +1 (white) for non-expressing cell types
  colors_1206$class_label <- makeQualitativePal(k_motifs, glasbey_use = T, skip = 0) # skip white color
  
  # Rename colors with pathway states -- including 0 for non-expressing cell types
  names(colors_1206$class_label) <- 0:(k_motifs -1)
  
  pathway_df$class_label = pathway_df$class_label  %>% as.character()
  # Make heatmap
  if(save_pdf){
    pheatmap(x_scaled,
             annotation_row = pathway_df %>% select(Cell_class, age, class_label, dataset),
             annotation_colors = colors_1206,
             show_colnames = F,
             show_rownames = F,
             clustering_method = clust_method,
             clustering_distance_rows = clust_dist_metric,
             treeheight_col = 0,
             cutree_rows = 12,
             fontsize = 10,
             color = magma(100),
             filename = paste(save_dir, '_global_dendrogram.pdf', sep = ""),
             height = 20,
             width = 4)
  }
  
  # Return the values used to compute global distance as well as pathway gene expression values and class labels
  return(list(x_scaled, pathway_df))
}
