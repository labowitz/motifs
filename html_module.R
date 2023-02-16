library(dplyr)
library(ggpubr)
library(ggplot2)
library(cluster)
library(Seurat)
library(pheatmap)
library(stylo)
library(tidyr)
library(ggsignif)
library(grid)
library(parallel)
library(gridExtra)
library(viridis)
library(superheat)
library(DIZutils)
library(stringr)
library(RColorBrewer)
library(zoo)
library(pracma)
library(DIZtools)

## Seurat object to import
master_seurat <- readRDS("./data/processed_data/master_seurat.RDS") # Most recent Seurat object
pca_proj <- readRDS("./data/processed_data/pca_proj.RDS") # PCA projection from 1928 HVGs
hvg_genes <- readRDS("./data/processed_data/hvg_genes.RDS") # list of HVGs (shouldn't be used)

## Set of files to import 
param_list <- readRDS("./data/processed_data/param_list.RDS") # List of pathway parameters
all_pathways <- readRDS("./data/processed_data/all_pathways.RDS") # List of all pathways

## Colors to import
colors_1206 <- readRDS("./data/processed_data/colors_1206.RDS") # List with labels and corresponding colors
glasbey <- readRDS("./data/processed_data/glasbey.RDS") # List of glasbey colors for categorical labels

black_pal <- function(n){
  colfunc <- colorRampPalette(c("white", "black"))
  return(colfunc(n))
}

makeQualitativePal <- function(n,               # Number of unique colors in palette
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
    col_vector = unlist(mapply(brewer.pal, 
                               qual_col_pals$maxcolors, 
                               rownames(qual_col_pals))
    ) # only 70 qualitative colors
    
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

## Read in master_seurat from csv
createSeurat <- function(file_path # path to the csv file 
){
  data <- read.csv(file_path, row.names=1)
  
  master_seurat <- CreateSeuratObject(counts=data)
  
  return(master_seurat)
}

# Makes a list out of the metadata_annotations stored in a csv
create_annotations <- function(file_path # csv with metadata_annotations
){
  annotations = read.csv(file_path)
  annotations_list = list()
  
  for (i in 1:ncol(annotations)){
    annotations_list[[i]]<-annotations[, i]
  }
  
  names(annotations_list) = colnames(annotations)
  
  return(annotations_list)
}

appendMetaData <- function(master_seurat, 
                           metadata_annotations
                           ){
  
  AddMetaData(master_seurat, metadata_annotations, col.names=NULL)
  
  ret
}

addGenes <- function(pathway_df = data.frame(),
                     pathway_name = "Notch",
                     pathway_genes = c()) {
  
  df <- data.frame(pathway = rep(pathway_name, length(pathway_genes)), gene = pathway_genes)
  pathway_df <- rbind(pathway_df, df)
  
  return(pathway_df)
}

getGenesList <- function(pathway_df = data.frame(),
                         pathway_name = "Notch"
                         ){
  if (!(pathway_name %in% pathway_df$pathway)) {
    stop("We don't know this pathway! Add the genes in this pathway to our list using the `addGenes` function!")
  }
  return(pathway_df %>% dplyr::filter(pathway == pathway_name) %>% pull(gene))
}

## Function to return the genes in a pathway stores in "all_pathways" object
genesPathway <- function(pathway_name = 'Notch',
                         pathway_df = data.frame(),
                         seurat_obj = master_seurat
){
  
  pathway_genes = getGenesList(pathway_name = pathway_name,
                              pathway_df = pathway_df)
  pathway_genes = pathway_genes[which(pathway_genes %in% row.names(seurat_obj))]
  
  return(pathway_genes)
}

## Retrieve the log-norm gene expression values and the annotations
makeMainDataFrame <- function(pathway_genes = c(), 
                              seurat_obj = master_seurat
){
  
  pathway_matrix <- FetchData(seurat_obj, 
                              pathway_genes) # Fetch log norm counts from Seurat object
  
  data_frame <- cbind(pathway_matrix, 
                      seurat_obj@meta.data %>%
                        dplyr::select(Tissue, 
                                      age, 
                                      dataset, 
                                      cell_ontology_class, 
                                      Cell_class, 
                                      cell_id) # Append metadata annotations
  )
  
  row.names(data_frame)<-1:dim(data_frame)[1] # Set the rownames
  
  data_frame$global_cluster = 1:dim(data_frame)[1] # here we haven't filtered anything. Comes directly from Seurat obj
  
  return(data_frame)
}

## MinMax Scaling of values
minMaxNorm <- function(x){
  maxs = apply(x, 2, max) # applies max operation to columns
  mins = apply(x, 2, min)
  
  # Scale the values
  for(i in 1:dim(x)[2]){
    x[, i] = (x[, i] - mins[i]) / (maxs[i] - mins[i])
  }
  
  return(x)
}

## Permute the rows for each column of data.frame
randomize_columns <- function(df = data.frame(),           # data.frame of gene expression values (rows are cells, columns are genes)
                              pathway_genes = c()  # Pathway genes list
){
  
  for(i in 1:length(pathway_genes))
    df[, pathway_genes[i]] = sample(df[, pathway_genes[i]])
  return(df)
  
}

## Returns MinMax normalized counts along with metadata at specified saturating quantile.
normalizedDevel <- function(seurat_obj = master_seurat,    # Which Seurat object to use
                            pathway_genes = c(),           # List of pathway genes
                            sat_val = 0.99,                # Saturating quantile
                            fill_zero_rows = F             # If a gene has all 0's, fill with very small number
){
  
  # Get the log norm gene counts and annotations
  data_frame <- makeMainDataFrame(pathway_genes, 
                                  seurat_obj)
  
  if(!'cell_id' %in% names(seurat_obj@meta.data))
    data_frame %>% mutate(cell_id = paste(global_cluster, 
                                          dataset, 
                                          sep = "_")
    ) -> data_frame
  
  # Get dataframe of only gene expression values in the subset of cell states we want
  x = data_frame[, pathway_genes]
  
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
  x <- minMaxNorm(x)
  
  # We can't do cosine distance on all zero genes so we just fill these with a very small value if we want
  if(fill_zero_rows)
    x[x == 0] = 10^-10
  
  # Modify the dataframe with the metadata to now have the MinMax scaled, saturated counts
  data_frame[, pathway_genes] <- x
  
  row.names(data_frame) <- data_frame$global_cluster # Reset column names
  
  return(data_frame)
}

## quickPipeline processes the counts from an object with a specified optimal number of clusters
quickPipeline <- function(seurat_obj = master_seurat,        # Seurat object
                          pathway_genes = c(),               # Pathway name to draw from all_pathways object
                          k_final = 25,                      # k value for clusters
                          min_genes_on = 1,                  # Min. number of genes for pathway to be "ON"
                          min_expr = 0.2,                    # Min. expression cutoff for gene to be "ON"
                          sat_val = 0.99
){
  
  # Get the MinMax scaled and normalized data and annotations
  data_frame <- normalizedDevel(seurat_obj = seurat_obj,
                                pathway_genes = pathway_genes,
                                sat_val = sat_val,
                                fill_zero_rows = F
  )
  
  # Compute the number of genes "ON" in the pathway
  data_frame$genes_on = rowSums(data_frame[, pathway_genes] > min_expr)
  
  # Filter out cell types with lower than desired expression
  data_frame %>% dplyr::filter(genes_on > min_genes_on) -> data_frame
  
  # Heatmap to for computing the distance tree for cell types
  p = pheatmap(data_frame[,pathway_genes],
               cutree_rows = k_final,
               clustering_method = 'ward.D2',
               clustering_distance_rows = dist.cosine(as.matrix(data_frame[, pathway_genes])),
               cluster_cols = F,
               silent = T,
               show_rownames = F,
               color = rev(colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name ="Greys")))(100))
  )
  
  # Get class labels
  cos_labels = cutree(p$tree_row, k = k_final) # Cosine tree at specified k-cut
  data_frame$class_label = cos_labels # Map class labels
  data_frame$class_label <- as.numeric(data_frame$class_label)
  
  data_frame %>% gather(key = 'gene', 
                        value ='expression', 
                        all_of(pathway_genes)) %>%
    group_by(class_label, 
             gene) %>%
    summarise(mean_expr = mean(expression), 
              n_cell_types = n()) %>%
    spread(gene, mean_expr) %>% 
    tibble::column_to_rownames(var = "class_label") -> x
  
  rownames(data_frame) <- rownames(data_frame)
  
  # Colors for pathway classes --- +1 (white) for non-expressing cell types
  colors_1206$class_label <- makeQualitativePal(n = k_final, 
                                                glasbey_use = T, 
                                                skip = 1) # skip white color
  
  names(colors_1206$class_label) <- data_frame$class_label %>% unique()
  
  p = pheatmap(data_frame[,pathway_genes],
               cutree_rows = k_final,
               annotation_row = data_frame %>% dplyr::select(class_label),
               annotation_colors = colors_1206,
               clustering_method = 'ward.D2',
               clustering_distance_rows = dist.cosine(as.matrix(data_frame[, pathway_genes])),
               cluster_cols = F,
               silent = T,
               show_rownames = F,
               color = rev(colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name ="Greys")))(100)),
               annotation_legend = F
  )
  
  # Return count matrix, annotations, and heatmap for selected pathway
  return(list('matrix' = x, 'data_frame' = data_frame, "plot" = p))
}

## Takes an array of different bootstrap runs and returns the average silhouette score for each k value from all runs
makeBootstrap_df <- function(s_boots){
  
  # rbind, so each column represents k value
  boot_mat  = do.call(rbind, s_boots)
  
  boot_df = data.frame(m = apply(boot_mat, 2, mean), s = apply(boot_mat, 2, sd)) # Mean over columns
  boot_df$k = 1:dim(boot_df)[1]
  
  return(boot_df)
}

## Function that computes silhouette scores for a given number of clusters
cluster_silhouette <- function(df = data.frame(), # data.frame in wide format with genes and class labels
                               pathway_genes = c(),# pathway genes
                               dist = 'euclidean',
                               return_singles = F # directly return the mean silhouette score
){
  
  # Assumes that the data.frame contains the pathway profiles in wide format
  # The input data.frame must contain global_cluster AND class_label
  ## class_label must be numeric friendly
  df %>% dplyr::select(c(all_of(pathway_genes), global_cluster, class_label)) -> df_mat
  
  # Silhouette requires numeric cluster labels
  labs <- df_mat$class_label %>% as.numeric()
  names(labs) <- df_mat$global_cluster
  x <- df_mat[, pathway_genes] %>% as.matrix
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

## Bootstraps silhouette scores from from global expression profiles
silhPathwayBootstrap <- function(pathway_genes = c(),                         # Pathway name
                                 seurat_obj = master_seurat,          # Seurat object
                                 k_max = 100,               # Maximum values of k_cutoffs to compute the silhouette score for
                                 clust_method = "ward.D2", # Clustering clust_method
                                 sat_val = 0.99,            # Saturating quantile
                                 dist_metric = 'cosine', # Distance metric
                                 pct_boots = 0.9,           # Fraction of cells we want to sample from dataset
                                 n_boots = 100,             # Number of bootstraps to run
                                 control_silh = F           # Whether this is a bootstrapping for clusters or a negative control that randomzes the data
){
  
  data_frame_main <- normalizedDevel(pathway_genes = pathway_genes,
                                     sat_val = sat_val,
                                     seurat_obj = seurat_obj, 
                                     fill_zero_rows = T
  )
  
  boot_list = list() # Bootstrap values list
  for(b in 1:n_boots){
    
    if(!control_silh){
      # Running bootstraps on the pathway genes
      # Filter out cell types that are not expressing.
      # but sample the list of cell types before to create a bootstrap distribution
      
      data_frame_main -> data_frame
      row.names(data_frame) <- data_frame$global_cluster
      
    }else{
      # Control bootstrap: re-shuffle the data to get a lower bound on silhouette score
      # by destroying gene-gene correlations but keeping the distribution of each gene.
      # Before, we still need to filter for the cell types to include

      # Randomize expression within those cells! We don't want gene expression values from cell types not originally included
      data_frame = randomize_columns(data_frame_main, pathway_genes)
      
    }
    
    x = data_frame[, pathway_genes] %>% as.matrix
    
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
      
      data_frame$class_label = clustering
      
      master_clustered = data_frame # data.frame with data and motif labels
      
      s1 <- cluster_silhouette(df = master_clustered,
                               pathway_genes = pathway_genes, # Pathway genes
                               dist = dist_metric # original clustering
      )
      
      # Compute silhouette score to return
      ss[k] = mean(s1$ms)
    }
    
    boot_list[[b]] = ss
  }
  
  return(boot_list)
}

## Compute pathway gene and randomized silhouette scores and plots them together
silhouettePlot <- function(pathway_genes = c(),       # Specify pathway genes
                           seurat_obj = master_seurat,
                           min_genes_on = 2,          # Min. number of genes for pathway to be "ON"
                           min_expr = 0.25,           # Min. expression cutoff for gene to be "ON"
                           n_bootstraps=10,           # Number of bootstrap replicates
                           pct_boots = 0.9,       # Percent of cells to ??
                           max_k = 100,               # Max. number of clusters to compute silhouette scores for
                           clust_metric = "cosine",   # Clustering metric
                           clust_method = "ward.D2"   # Clustering method
){
  
  # 1. Compute the real pathway silhouette scores
  s_boots = silhPathwayBootstrap(pathway_genes = pathway_genes,
                                 seurat_obj = seurat_obj,
                                 k_max= max_k,
                                 dist_metric = clust_metric,
                                 clust_method = clust_method,
                                 control_silh = F,
                                 n_boots = n_bootstraps,
                                 pct_boots = pct_boots
  )
  
  boot_df = makeBootstrap_df(s_boots)
  
  # 2. Compute the negative control silhouette scores
  s_boots = silhPathwayBootstrap(pathway_genes = pathway_genes,
                                 seurat_obj = seurat_obj,
                                 k_max= max_k,
                                 dist_metric = clust_metric,
                                 clust_method = clust_method,
                                 control_silh = T,
                                 n_boots = n_bootstraps,
                                 pct_boots = pct_boots
  )
  
  boot_df_control = makeBootstrap_df(s_boots)
  
  # 3. Plot together
  boot_df %>% ggplot(aes(x = k, y = m)) +
    geom_ribbon(aes(ymin = m - s, 
                    ymax = m + s), 
                fill = 'lightblue', 
                alpha = 0.2) +
    geom_line(color ='blue') + 
    geom_ribbon(data = boot_df_control, 
                aes(ymin = m - s, 
                    ymax = m + s), 
                alpha = 0.2) +
    geom_line(data = boot_df_control) + 
    theme_pubr(base_size = 16) + 
    ylab('Silhouette') +
    xlab('Number clusters') -> g
  
  # Returns plot and both data frames
  return(list(g, boot_df, boot_df_control, s_boots))
}

## Computes and plots Z-score based on bootstrapped pathway and randomized expression values
silhouette_zscore <- function(silh_result,                        # Result from running silhouettePlot function
                              min_expr = 0.2,                     # Min. expression cutoff for gene to be "ON"
                              x_offset = 1,                       # X-axis of k (number of clusters) begins plotting at x_offset + 1
                              min.y = 0.1,                        # Min. cut-off on y-axis of the silhouette score plot
                              max.y = 0.7,                         # Max. cut-off on y-axis of the silhouette score plot
                              k_max = 100
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
    geom_vline(xintercept = which.max(z_score), 
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
  
  return(list(g, which.max(z_score)))
}

## Rename
## Compute silhouette score for given labels, and returns mean silhouette score per class
clusteringSilhouette <- function(df = data.frame(),         # data.frame in wide format with genes and class labels
                                 dist_metric = 'cosine',    # distance metric
                                 pathway_genes = c()        # List of genes in pathway
){
  # Compute the silhouette score for k cluster labels
  if(dist_metric == 'cosine'){
    d = silhouette(df$class_label %>% as.numeric, 
                   dist = dist.cosine(df[, pathway_genes] %>% as.matrix)
    )
    
  }else{
    d = silhouette(df$class_label %>% as.numeric, 
                   dist = dist(df[, pathway_genes] %>% as.matrix)
    )
    
  }
  
  # Store silhouette results in data.frame of class label, neigbhoring cluster, and silhouette score
  silh_df = data.frame(class_label = d[, 1], 
                       neighbor = d[, 2], 
                       silh = d[, 3]
  )
  
  # Take the silhouette score for a cluster as the average of that of all its members
  silh_df %>% 
    group_by(class_label) %>% 
    summarise(mean_silh = mean(silh)) -> silh_df
  
  return(silh_df)
}

## Computes pathway and global diversity, along with mean ss for each class
globalClustering <- function(df_devel = data.frame(),   # data.frame in wide format with genes and class labels
                             pathway_genes = c(),       # List of pathway genes
                             seurat_obj = c(),          # Seurat object
                             k_final = 25,              # k value for clusters
                             n_pcs = 100,               # Number of PCs to use
                             manual_embedding = c(),    # Matrix should represent a latent space (i.e. PCA, UMAP), do not provide if using n_pcs
                             dist_metric = 'euclidean'  # Distance metric
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
    umap_coords <- Embeddings(seurat_obj, reduction = 'pca')
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
                      annotation_col = df_devel %>% 
                        dplyr::select(age, dataset, Tissue, Cell_class, class_label),
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
      return(dist_class[upper.tri(dist_class)] %>% mean())
      
    }else{      # If there is only one member in the class, then pairwise similarity is zero
      return(0)
    }
    
  }) -> diversity_pathway
  
  # Perform the same pairwise similarity computation on global labels
  
  global_labels = cutree(p_global$tree_col, 
                         k = k_final)
  df_devel$global_label = global_labels
  
  sapply(1:k_final, function(x){
    class_cells = df_devel$cell_id[df_devel$global_label == x]
    
    if(length(class_cells)>1){
      dist_mat[class_cells, class_cells] -> dist_class
      return(dist_class[upper.tri(dist_class)] %>% mean())
      
    }else{
      return(0)
    }
    
  }) -> diversity_global
  
  # Transcriptome diversity
  df_devel %>% 
    group_by(global_label) %>% 
    count %>% 
    as.data.frame() -> global_stats
  global_stats$diversity = diversity_global
  
  # Pathway diversity
  df_devel %>% 
    group_by(class_label) %>% 
    count %>% 
    as.data.frame() %>%
    mutate(class_label = as.numeric(class_label)) %>%
    arrange(class_label) -> path_stats
  path_stats$diversity = diversity_pathway
  
  # Compute silhouette scores for this clustering
  # By default, we perform gene clustering using the cosine distance and ward.D2
  # so we also use cosine distance for the silhouette score
  
  pathway_silh = clusteringSilhouette(df = df_devel, 
                                      dist_metric = 'cosine', 
                                      pathway_genes = pathway_genes)
  
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
positiveControl <- function(df_devel = data.frame(),
                            seurat_obj = c(),
                            n_pcs = 100,
                            n_random = 100,
                            manual_embedding = c(),
                            dist_metric  = 'euclidean'
){
  # Cells are filtered
  df_devel$class_label <- df_devel$class_label %>% as.character()
  k_final = length(df_devel$class_label %>% unique)
  
  # PCA embedding
  # Number of principal components
  # umap_coords by default uses cell_id as the row.name
  if(length(manual_embedding) == 0){
    umap_coords <- Embeddings(seurat_obj,
                              reduction = 'pca')
    umap_coords <- umap_coords[, 1:n_pcs]
    
  }else{
    # Here the latent space already has the row.names
    # dimensions should be specified before
    umap_coords = manual_embedding
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
        
        return(dist_class[upper.tri(dist_class)] %>% mean())
        
      }else{
        return(0)
      }
      
    }) -> diversity_pathway
    
    # Pathway
    df_devel %>% 
      group_by(class_label) %>% 
      count %>% 
      as.data.frame() %>%
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

parallelPipeline <- function(pathway_genes=c(),          # Pathway list
                             k_final = c(),              # k value for clusters
                             seurat_obj = c(),           # Seurat object
                             manual_cell_types = c(),    # Any cell types to filter out
                             n_pcs = 100,                # Number of PCs to use if using PC embedding
                             manual_embedding = c(),       # Else user-specified Embedding matrix
                             dist_metric = 'euclidean',  # Distance metric to use to compute similarities
                             min_genes_on=1,
                             min_expr=0.2
){
  
  # quickPipeline takes manual_cell_types as the list with cell_ids to include
  # If empty, the function will filter cell_types based on the expression of the selected genes
  res_list = quickPipeline(pathway_genes = pathway_genes,
                           seurat_obj = seurat_obj,
                           k_final = k_final,
                           min_genes_on=min_genes_on,
                           min_expr=min_expr
  )
  
  stats_list = globalClustering(df_devel = res_list$data_frame,
                                seurat_obj = seurat_obj,
                                k_final = k_final,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric,
                                pathway_genes = pathway_genes
  )
  return(stats_list)
}

# Runs negative control tests for cell-type specific profiles by computing
# __ and diversity from non-pathway genes
controlDiversity <- function(null_dist = c(),
                             pathway_genes = c(),
                             k_final = c(),
                             seurat_obj= c() ,
                             n_rand=100,
                             filter_manual_cells = c(),
                             verbose = F, 
                             n_pcs = 100,
                             manual_embedding = c(),
                             dist_metric = 'euclidean',
                             min_genes_on=1,
                             min_expr=0.2
){
  
  control_stats = list()
  
  rand_list = list()
  for(rr in 1:n_rand)
    rand_list[[rr]] = sample(null_dist, length(pathway_genes))
  
  # Run control pipeline in parallel for the whole list of random pathways
  # specify the manual list of cell types: overrides the filtering by min.expression
  control_stats = mclapply(rand_list, 
                           parallelPipeline,
                           k_final = k_final,
                           seurat_obj = seurat_obj,
                           manual_cell_types = filter_manual_cells,
                           n_pcs = n_pcs,
                           manual_embedding = manual_embedding,
                           dist_metric = dist_metric,
                           min_genes_on=min_genes_on,
                           min_expr=min_expr,
                           mc.cores = 10)
  
  # make control_df
  print(control_stats[[1]])
  do.call(rbind, lapply(control_stats, 
                        function(x){rbind(x$pathway %>% 
                                            dplyr::select(-rank), 
                                          x$global)})) -> control_df
  
  return(control_df)
  
}


fullControlPathway <- function(pathway_genes = c(),      # List of pathway genes
                               k_final = c(),            # k value for clusters
                               seurat_obj = c(),         # Seurat object
                               null_list = c(),          # List of genes for control, HVGs to generate random data
                               n_samples = 100,
                               filter_manual = F,
                               min_expr = 0.2,           # Min. expression cutoff for gene to be "ON"
                               min_genes_on = 1,         # Min. number of genes for pathway to be "ON"
                               n_pcs = 100,              # Number of PCs to use
                               manual_embedding = c(),   # Number of principal components to consider for diversity metric
                               dist_metric = 'euclidean' # distance metric to be use in the global embedding (e.g., PCA) for global similarity across cell types
){
  
  # 1. Run the quick pipeline for the real pathway -- does not use any global metric for transcriptome
  #    with the specified k (found previously as optimal)
  res_list = quickPipeline(pathway_genes = pathway_genes,
                           seurat_obj = seurat_obj,
                           k_final = k_final, 
                           min_genes_on = min_genes_on, 
                           min_expr = min_expr
  )
  
  # 2. Compute diversity for the cell types containing each of the pathway profiles
  # Here we use an embedding to compute the global metric (PCA, VAE, ICA, etc)
  stats_list = globalClustering(df_devel = res_list$data_frame,
                                seurat_obj = seurat_obj,
                                k_final = k_final,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric,
                                pathway_genes = pathway_genes)
  
  
  # 3. Negative control: random sets of genes from null list
  # Note: here we want two types of filtering:
  #       current: filter cell types based on the expression of the random set
  #       probably better: filter cell types that express the real pathway so the diversity is comparable
  # new option to filter the same cell types as the real pathway
  if(filter_manual){
    manual_cell_types = res_list$data_frame$cell_id
  }else{
    manual_cell_types = c()
  }
  
  # Add the filter using filter_manual_cells parameter
  control_df = controlDiversity(null_dist = null_list,
                                pathway_genes = pathway_genes,
                                k_final = k_final,
                                seurat_obj = seurat_obj,
                                n_rand = n_samples,
                                filter_manual_cells = manual_cell_types,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric,
                                min_expr = min_expr,     # Min. expression cutoff for gene to be "ON"
                                min_genes_on = min_genes_on)
  
  # Split negative control results: pathway (random genes), transcriptome (global diversity)
  control_df_path = control_df %>% dplyr::filter(type == 'pathway')
  control_df_global = control_df %>% dplyr::filter(type == 'transcriptome')
  
  # 4. Positive control: re-distribute pathway profiles in randomly selected cell types.
  # returns diversity scores only for the fake pathway (no global calculation)
  pos_control = positiveControl(df_devel=res_list$data_frame,
                                seurat_obj = seurat_obj,
                                n_random = n_samples,
                                n_pcs = n_pcs,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric
  )
  # Merge all in data.frame
  pos_control = do.call(rbind, pos_control)
  
  
  df_diversity = rbind(data.frame(d = control_df_path$diversity, 
                                  type = 'null dist'), #random sets
                       data.frame(d = control_df_global$diversity, 
                                  type = 'transcriptome dist'), #global dendrogram using same cells from random sets
                       data.frame(d = stats_list$pathway$diversity, 
                                  type = 'pathway'), # actual pathway
                       data.frame(d = stats_list$global$diversity, 
                                  type = 'transcriptome'), # global dendrogram using the same cells than the real pathway
                       data.frame(d = pos_control$diversity, 
                                  type = 'pos control')) #randomized pathway profiles across cell types
  
  df_recurrence = rbind(data.frame(d = control_df_path$diversity * control_df_path$n, 
                                   type = 'null dist'),
                        data.frame(d = control_df_global$diversity * control_df_global$n, 
                                   type = 'transcriptome dist'),
                        data.frame(d = stats_list$pathway$diversity * stats_list$pathway$n, 
                                   type = 'pathway'),
                        data.frame(d = stats_list$global$diversity * stats_list$global$n, 
                                   type = 'transcriptome'),
                        data.frame(d = pos_control$diversity * pos_control$n, 
                                   type ='pos control'))
  
  # We can also return the main data.frame, with class_label, diversity and rank for each profile
  clustered_data <- res_list$data_frame %>% left_join(stats_list$pathway %>%
                                                        rename(class_label = label) %>%
                                                        dplyr::select(rank, 
                                                                      diversity, 
                                                                      class_label, 
                                                                      n), 
                                                      by = 'class_label')
  
  return(list(diversity = df_diversity, 
              recurrence = df_recurrence, 
              profiles = clustered_data, 
              rank = stats_list)
  )
}

## ECDF plotting function
ecdfDiversity <- function(control_res = list(),
                          plot_null = F,
                          pathway_name = 'Bmp'
){
  
  if(!plot_null){
    df = control_res$diversity %>% 
      dplyr::filter(type %in% c('pathway','transcriptome','pos control'))
    
    lines_palette = c("#00AFBB", "#000000", "#a8a8a8") # color order: pathway, transcriptome, random
    
  }else{
    df = control_res$diversity %>% dplyr::filter(type %in% c('pathway','transcriptome','pos control','null dist'))
    lines_palette = c("#00AFBB", "#000000", "#a8a8a8",'#6f7e96') # color order: pathway, transcriptome, random
    
  }
  quantile_line = quantile(df %>% 
                             dplyr::filter(type=='transcriptome') %>% 
                             pull(d), 
                           0.9)
  
  df  %>% ggplot(aes(x = d, col = type)) + 
    stat_ecdf()  + 
    ggtitle("ECDF of Pathway Dispersion") +
    scale_x_continuous(trans = 'sqrt') + 
    theme_pubr(base_size = 16) +
    scale_color_manual(values = lines_palette) +
    ylab('Fraction of clusters') + 
    xlab('Cell state dispersion') +
    geom_hline(yintercept = 0.95,
               linetype="dashed", 
               color = "lightgray", 
               linewidth = 2) -> ecdf_plot
  
  return(ecdf_plot)
}

## Plotting function for ranked diversity values
diversityPlot <- function(stats_list = list(), # Pathway and global diversity statistics from running fullControlPathway
                          title_main = 'Notch' # Pathway name for plot title
){
  # Extract pathway and global diversity values sorted by rank
  path_stats = stats_list$pathway
  global_stats = stats_list$global
  upper_quantile = 0.75
  
  g <- path_stats %>% ggplot(aes(x = rank, y = diversity)) + 
    geom_point() +
    geom_hline(yintercept = quantile(global_stats$diversity, 
                                     upper_quantile), 
               linetype = "dashed", 
               color = "lightgray", 
               linewidth = 1) +
    geom_hline(yintercept = quantile(global_stats$diversity, 
                                     1 - upper_quantile), 
               linetype = "dashed", 
               color = "lightgray", 
               linewidth = 1) +
    ylab("Cell type diversity (PCA)") + 
    theme_pubr(base_size = 16) +
    xlab("Pathway profile") + 
    ggtitle(title_main) +
    annotate(geom = "text", 
             x = 4, 
             y = quantile(global_stats$diversity, 
                          upper_quantile) + 2, 
             label = "Expected",
             size=16/.pt,
             color="black")
  
  return(g)
}## Plotting function for ranked diversity values
diversityPlot <- function(stats_list = list(), # Pathway and global diversity statistics from running fullControlPathway
                          title_main = 'Notch' # Pathway name for plot title
){
  # Extract pathway and global diversity values sorted by rank
  path_stats = stats_list$pathway
  global_stats = stats_list$global
  upper_quantile = 0.75
  
  g <- path_stats %>% ggplot(aes(x = rank, y = diversity)) + 
    geom_point() +
    geom_hline(yintercept = quantile(global_stats$diversity, 
                                     upper_quantile), 
               linetype = "dashed", 
               color = "lightgray", 
               linewidth = 1) +
    geom_hline(yintercept = quantile(global_stats$diversity, 
                                     1 - upper_quantile), 
               linetype = "dashed", 
               color = "lightgray", 
               linewidth = 1) +
    ylab("Cell type diversity (PCA)") + 
    theme_pubr(base_size = 16) +
    xlab("Pathway profile") + 
    ggtitle(title_main) +
    annotate(geom = "text", 
             x = 4, 
             y = quantile(global_stats$diversity, 
                          upper_quantile) + 2, 
             label = "Expected",
             size=16/.pt,
             color="black")
  
  return(g)
}

motif_heatmap <- function(control_res=data.frame(),
                          pathway_genes=c(),
                          diverse_quantile=0.90
){
  ## Plot the most diverse profiles, aka the motifs
  # 1. Select all profiles 
  diverse_df <- control_res$profiles
  # 2. tidy data.frame to average gene expression within a pathway profile 
  diverse_df %>% pivot_longer(cols = pathway_genes, 
                              names_to = 'gene', 
                              values_to = 'expression') %>% 
    select(cell_id, cell_ontology_class, Tissue, Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy 
  
  diverse_tidy %>% group_by(class_label,gene) %>% summarise(mean_expr = mean(expression), 
                                                            rank = mean(rank),diversity = mean(diversity), 
                                                            cell_types = mean(n)) -> diverse_tidy
  
  # 3. wide format 
  diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types), names_from=gene,values_from=mean_expr) %>% tibble::column_to_rownames('class_label')-> diverse_mat
  
  # 4. We do the filtering here either for motifs or for non-diverse profiles 
  control_res$diversity %>% 
    dplyr::filter(type=='transcriptome') %>% # choose the null model 
    pull(d) %>% quantile(0.90) -> diverse_quantile
  
  diverse_mat %>% dplyr::filter(diversity>diverse_quantile) -> diverse_mat 
  
  motif_heatmap <- superheat(diverse_mat[,pathway_genes],
                             pretty.order.rows = T,
                             heat.pal = black_pal(10),
                             bottom.label.text.angle = 90, 
                             yr = sqrt(diverse_mat$cell_types),
                             yr.plot.type='bar',
                             yr.axis.name = "N cell types",
                             row.title = "Pathway motifs",
                             column.title = "Pathway genes",
                             bottom.label.size = 0.3,
                             grid.hline.col = "white",
                             grid.hline.size = 2, 
                             
  )
  
  return(motif_heatmap)
}

motif_ct_heatmap <- function(control_res=data.frame(),
                             pathway_genes=c(),
                             diverse_quantile=0.90,
                             type="motif"
){
  
  div_res <- diverseFilt(control_res = control_res,
                         pathway_genes = pathway_genes,
                         diverse_quantile = diverse_quantile,
                         type=type
  )
  
  # 6. Tissue and organ distributions of selected profiles
  div_res$diverse_df %>% 
    dplyr::filter(class_label %in% row.names(div_res$diverse_mat) ) %>% 
    dplyr::select(Tissue, class_label) %>% 
    dplyr::filter(!Tissue %in% c('mat','endoderm')) %>% 
    group_by(class_label,Tissue) %>% 
    count  %>% 
    pivot_wider(id_cols=class_label,
                names_from = Tissue,
                values_from = n,
                values_fill = 0) -> tissue_distribution
  
  
  x = tissue_distribution %>% 
    ungroup %>% 
    dplyr::select(-class_label) %>% 
    as.matrix() 
  row.names(x) <- tissue_distribution$class_label
  # make tissue names pretty 
  colnames(x)<-str_replace(string = colnames(x),pattern ="_",replacement = ' ') %>% firstup() 
  
  x <- x[,sort(colnames(x))]
  
  data.frame(x) %>% 
    rownames_to_column("profiles") %>% 
    gather("tissues", values, -profiles) -> x
  
  x$values <- sqrt(x$values)
  
  plt <- ggplot(x, 
                aes(tissues,profiles,fill=values)) + 
    geom_tile(color="grey") + 
    scale_y_discrete(position = "right", limits=rev) + 
    scale_fill_distiller(palette = "PuRd", direction=1) + 
    guides(fill=guide_colorbar("sqrt(No.\ncell types)")) +
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank()) + 
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_text(angle=45, 
                                     vjust = 1.1, 
                                     hjust = 1),
          axis.title=element_blank(),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10),
          plot.title = element_text(size=20))
  return(plt)
}

## Computes the global transcriptome distance using a specified embedding and plots pathway class labels on global dendrogram
global_dendr <- function(control_res = list(),                         # Main results object
                         seurat_obj = master_seurat,                             # Seurat object
                         hvg_genes = c(),                              # List of genes to create the global dendrogram
                         dist_metric ='cosine',                        # Clustering distance in global space. applies for expression and PCA spaces
                         clust_method = "ward.d2",                     # Hierarchical clustering method
                         use_pca = F,                                  # Whether to use PCA coordinates instead of gene expression for global dendrogram
                         n_pcs = 1:100                                 # If use_pca, then n_pcs to use
){
  
  # Get the pathway states of the cells which have the pathway "ON"
  profiles_df = control_res$profiles
  row.names(profiles_df) <- profiles_df$cell_id
  
  meta_full <- master_seurat@meta.data %>% 
    dplyr::select(cell_ontology_class,
                  Cell_class, 
                  dataset, 
                  age)
  
  meta_full$cell_id <- row.names(meta_full)
  
  # Merge the pathway states with the whole dataset
  meta_full$class_label = 0 # assign cells with pathway "OFF" a class_label = 0
  meta_full[row.names(profiles_df), ]$class_label <- profiles_df$class_label
  
  profiles_df <- meta_full
  
  # Compute the global distance in the specified Embedding space
  if(!use_pca){
    # Use highly variable genes
    hvg_df = makeMainDataFrame(hvg_genes, 
                               master_seurat = seurat_obj)
    hvg_genes = hvg_genes[hvg_genes %in% colnames(hvg_df)]
    
    # Filter only for expressing cell types in this pathway
    hvg_df %>% dplyr::filter(cell_id %in% profiles_df$cell_id) -> hvg_df
    row.names(hvg_df) <- hvg_df$cell_id
    
    # Scale the gene expression values so mean is zero
    x = hvg_df[, hvg_genes]
    x_scaled = scale(x) # standardScaler
    row.names(x_scaled) <- hvg_df$cell_id
    
  }else{
    # Get the values in PCA space
    x_scaled = Embeddings(seurat_obj, reduction='pca')
    x_scaled = x_scaled[profiles_df$cell_id, n_pcs] # already set row names as cell id
    
  }
  
  # Compute distance in the global transcriptome space
  if(dist_metric == 'cosine') clust_dist_metric = dist.cosine(x_scaled) else clust_dist_metric = dist(x_scaled)
  
  k_motifs = length(profiles_df$class_label %>% unique) # number of pathway states
  
  # Colors for pathway classes --- +1 (white) for non-expressing cell types
  colors_1206$class_label <- makeQualitativePal(k_motifs, 
                                                glasbey_use = T, 
                                                skip = 1) # skip white color
  
  # Rename colors with pathway states -- including 0 for non-expressing cell types
  names(colors_1206$class_label) <- 0:(k_motifs -1)
  
  profiles_df$class_label = profiles_df$class_label %>% as.character()
  # Make heatmap
  heatmap_plot <- pheatmap(x_scaled,
                           annotation_row = profiles_df %>% 
                             dplyr::select(Cell_class, 
                                           age, 
                                           class_label, 
                                           dataset),
                           annotation_colors = colors_1206,
                           show_colnames = F,
                           show_rownames = F,
                           clustering_method = clust_method,
                           clustering_distance_rows = clust_dist_metric,
                           treeheight_col = 0,
                           cutree_rows = 12,
                           fontsize = 10,
                           color = magma(100)
  )
  
  # Return the values used to compute global distance as well as pathway gene expression values and class labels
  return(list(plt = heatmap_plot, profiles_df = profiles_df))
}

global_umap <- function(control_res = list(),                         # Main results object
                        seurat_obj = master_seurat,                   # Seurat object
                        hvg_genes = c(),                              # List of genes to create the global dendrogram
                        dist_metric ='cosine',                        # Clustering distance in global space. applies for expression and PCA spaces
                        clust_method = "ward.d2",                     # Hierarchical clustering method
                        use_pca = F,                                  # Whether to use PCA coordinates instead of gene expression for global dendrogram
                        n_pcs = 1:100
){
  # Get the pathway states of the cells which have the pathway "ON"
  profiles_df = control_res$profiles
  row.names(profiles_df) <- profiles_df$cell_id
  
  meta_full <- master_seurat@meta.data %>% 
    dplyr::select(cell_ontology_class,
                  Cell_class, 
                  dataset,
                  age)
  
  meta_full$cell_id <- row.names(meta_full)
  
  # Merge the pathway states with the whole dataset
  meta_full$class_label = 0 # assign cells with pathway "OFF" a class_label = 0
  meta_full[row.names(profiles_df), ]$class_label <- profiles_df$class_label
  
  profiles_df <- meta_full
  
  # Compute the global distance in the specified Embedding space
  if(!use_pca){
    # Use highly variable genes
    hvg_df = makeMainDataFrame(hvg_genes, 
                               master_seurat = seurat_obj)
    hvg_genes = hvg_genes[hvg_genes %in% colnames(hvg_df)]
    
    # Filter only for expressing cell types in this pathway
    hvg_df %>% dplyr::filter(cell_id %in% 
                               profiles_df$cell_id) -> hvg_df
    row.names(hvg_df) <- hvg_df$cell_id
    
    # Scale the gene expression values so mean is zero
    x = hvg_df[, hvg_genes]
    x_scaled = scale(x) # standardScaler
    row.names(x_scaled) <- hvg_df$cell_id
    
  }else{
    # Get the values in PCA space
    x_scaled = Embeddings(seurat_obj, reduction='pca')
    x_scaled = x_scaled[profiles_df$cell_id, n_pcs] # already set row names as cell id
    
  }
  
  # Compute distance in the global transcriptome space
  if(dist_metric == 'cosine') clust_dist_metric = dist.cosine(x_scaled) else clust_dist_metric = dist(x_scaled)
  
  k_motifs = length(profiles_df$class_label %>% unique) # number of pathway states
  
  # Colors for pathway classes --- +1 (white) for non-expressing cell types
  colors_1206$class_label <- makeQualitativePal(k_motifs, 
                                                glasbey_use = T, 
                                                skip = 0) # skip white color
  
  # Rename colors with pathway states -- including 0 for non-expressing cell types
  names(colors_1206$class_label) <- 0:(k_motifs -1)
  
  profiles_df$class_label = profiles_df$class_label  %>% as.character()
  
  plt <- DimPlot(seurat_obj, 
                 reduction = "umap", 
                 cells.highlight = profiles_df %>% 
                   filter(class_label != "0") %>%
                   rownames(), 
                 cols.highlight = alpha("black", 0.5), 
                 col = alpha("grey", 0.5)) + 
    theme_void() + 
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", 
                                      fill=NA, 
                                      linewidth = 0.5)) +
    ggtitle("Global UMAP: Pathway ON")
  
  return(plt)
}

diversity_plot <- function(stats_list = list() # Pathway and global diversity statistics from running fullControlPathway
){
  # Extract pathway and global diversity values sorted by rank
  path_stats = stats_list$pathway
  global_stats = stats_list$global
  upper_quantile = 0.75
  
  g <- path_stats %>% 
    ggplot(aes(x = rank, y = diversity)) + 
    geom_point() +
    geom_hline(yintercept = quantile(global_stats$diversity, upper_quantile), 
               linetype = "dashed", 
               color = "lightgray", 
               size = 1) +
    geom_hline(yintercept = quantile(global_stats$diversity, 1 - upper_quantile), 
               linetype = "dashed", 
               color = "lightgray", 
               size = 1) +
    ylab("Cell type diversity (PCA)") + 
    theme_pubr(base_size = 10) +
    xlab("Pathway profile") + 
    ggtitle('Profile Dispersion Plot') +
    annotate(geom = "text", 
             x = 3, 
             y = quantile(global_stats$diversity, 
                          upper_quantile) + 2, 
             label = "Expected", 
             color="black")
  
  return(g)
}

## Actually compute the diversity values of the specified number of profiles and plot them
rank_diversity <- function(pathway_genes =c(),
                           k_final = 20,
                           min_expr = 0.2,
                           min_genes_on = 2,
                           manual_embedding = c(),
                           dist_metric = 'euclidean',
                           make_plot = T,
                           seurat_obj=master_seurat
){
  
  # takes the pathway name as input:
  # all parameters come from upstream functions
  res_list = quickPipeline(pathway_genes = pathway_genes,
                           seurat_obj = seurat_obj,
                           k_final = k_final, 
                           min_genes_on = min_genes_on, 
                           min_expr = min_expr)
  
  # Get the diversity values
  stats_list = globalClustering(df_devel = res_list$data_frame,
                                seurat_obj = seurat_obj,
                                k_final = k_final,
                                manual_embedding = manual_embedding,
                                dist_metric = dist_metric,
                                pathway_genes = pathway_genes)
  
  if(make_plot){
    return(diversity_plot(stats_list = stats_list))
  }else{
    return(stats_list)
  }
}


