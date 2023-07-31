source("./scripts/analysis/imports.R")
library(MASS)

# Function to compute the cosine distance of two gene expression profiles
computePathwayDist <- function(pipe_run){
  
  pathway_dist = dist.cosine(as.matrix(pipe_run$data_frame[,pathway_genes]))
  
  pathway_dist_df <- data.frame(t(combn(pipe_run$data_frame$cell_id,2)), 
                                as.numeric(pathway_dist))
  
  colnames(pathway_dist_df) <- c("cell_id_1", "cell_id_2", "pathway_dist")
  
  return(pathway_dist_df)
  
}

# Computes the pathway, global, and null dist. of profile distances in expressing cell types
dfPathwayGlobalDist <- function(pathway_genes = c(),
                                  pathway_df = pathway_df,
                                  seurat_obj = master_seurat,
                                  k_final = optimal_k_pathway,
                                  min_genes_on = min_genes_pathway,
                                  min_expr = min_expr_threshold){
  
  pipe_run <- quickPipeline(seurat_obj = seurat_obj,
                            pathway_genes = pathway_genes,
                            k_final = optimal_k_pathway, 
                            min_genes_on = min_genes_on, 
                            min_expr = min_expr
  )
  
  scrambled_pipe_run <- pipe_run
  
  scrambled_pipe_run$data_frame[,pathway_genes] <- randomizeColumns(df = pipe_run$data_frame[,pathway_genes],
                                                                    pathway_genes = pathway_genes)
  
  real_pathway_dist_df <- computePathwayDist(pipe_run)
  
  scrambled_pathway_dist_df <- computePathwayDist(scrambled_pipe_run)
  colnames(scrambled_pathway_dist_df)[colnames(scrambled_pathway_dist_df) == 'pathway_dist'] <- 'random_dist'
  
  # GLOBAL DISTANCE: euclidean or cosine distance in PCA space
  global_coords = Embeddings(seurat_obj, reduction='pca')
  global_coords = global_coords[pipe_run$data_frame$cell_id, n_pcs] # already set row names as cell id
  
  global_dist = dist(global_coords)
  
  global_dist_df <- data.frame(t(combn(rownames(global_coords),2)), 
                               as.numeric(global_dist))
  
  colnames(global_dist_df) <- c("cell_id_1", "cell_id_2", "global_dist")

  df <- merge(x = real_pathway_dist_df, 
              y = scrambled_pathway_dist_df,
              by.x = c("cell_id_1", 
                       "cell_id_2"))
  
  df <- merge(x = df, 
              y = global_dist_df, 
              by.x = c("cell_id_1", 
                       "cell_id_2"))
  colnames(df) <- c("cell_id_1", "cell_id_2", "pathway_dist", "random_dist", "global_dist")
  
  df <- merge(df, 
              pipe_run$data_frame %>% 
                dplyr::select(c(cell_id, class_label)), 
                by.x = "cell_id_1", 
                by.y = "cell_id")
  
  colnames(df) <- c("cell_id_1", "cell_id_2", "pathway_dist", "random_dist", "global_dist", "class_label_id_1")
  
  df <- merge(df, 
              pipe_run$data_frame %>% 
                dplyr::select(c(cell_id, class_label)), 
              by.x = "cell_id_2", 
              by.y = "cell_id")
  
  colnames(df) <- c("cell_id_1", "cell_id_2", "pathway_dist", "random_dist", "global_dist", "class_label_id_1", "class_label_id_2")
  
  return(df)
}

# Ran the silhouette score and dispersion computations earlier and stored here
silh_res_dir = "./scripts/figures/peak_analysis/silhouette_res/silh_rds/alejo_res/"
dispersion_dir = "./scripts/figures/peak_analysis/dispersion/disp_rds/"

# Import data -- lists of genes 
pathway_df <- read.table("./data/raw_data/pathbank/pathway_df_linux_format.tsv", 
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
# fix the names of some pathways for Linux formatting
pathway_df$pathway <- pathway_df$pathway %>% str_replace('/', ' ')
pathway_df$pathway <- pathway_df$pathway %>% str_replace('\\(', ' ')
pathway_df$pathway <- pathway_df$pathway %>% str_replace('\\)', ' ')
pathway_df$pathway <- pathway_df$pathway %>% str_replace('–', '-')

fig_5_df <- read.csv("scripts/figures/peak_analysis/silhouette_res/dispersion_states_figure5.csv")

pathway_df <- pathway_df[pathway_df$pathway %in% fig_5_df$pathway_name,]

## Save pathway data.frame
saved_files <- paste(list.files(path=silh_res_dir,
                                pattern = "_ale.RDS"), 
                     sep ="") # ale is Alejandro's personal ID 

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
pathway_list <- full_pathway_list[saved_idx]

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.3 # tgfb minimum expression threshold for gene to be on
optimal_k_pathway = 30 # doesn't matter for this!
diverse_quantile = 0.9

n_pcs = 1:30

dataset_list <- c("1m 10x",
                  "3m 10x",
                  "18m 10x",
                  "21m 10x",
                  "24m 10x",
                  "30m 10x")

for (i in pathway_list){
  
  if(i %in% c('Notch receptors, Dll ligands and Fringe proteins',
              'Tgf-beta family receptors')){
    min_expr_threshold = 0.2
  }
  else{
    min_expr_threshold = 0.3
  }
  
  pathway_genes = genesPathway(pathway_name = i,
                               pathway_df = pathway_df,
                               seurat_obj = master_seurat)
  
  
  seurat_obj = subset(x = master_seurat, subset = (dataset %in% dataset_list))
  
  dfPathwayGlobalDist(pathway_genes = pathway_genes, 
                      seurat_obj = seurat_obj) -> df
  
  df %>% 
    dplyr::filter(global_dist >40 & global_dist < 90) %>%  
    ggplot(aes(x = pathway_dist)) + 
    geom_density() + 
    geom_density(aes(x = random_dist), color = 'red')
  
  # non disperse
  df %>% 
    dplyr::filter(global_dist <40 ) %>%  
    ggplot(aes(x = pathway_dist)) + 
    geom_density() + 
    geom_density(aes(x = random_dist), color = 'red')
  
  df %>% 
    pivot_longer(cols = c("pathway_dist", "random_dist"), 
                 names_to = "type", 
                 values_to = "distance") -> adult_tidy
  
  adult_tidy %>% 
    dplyr::filter(global_dist >40 & global_dist <90 ) %>% 
    ggplot(aes(x = distance, color = type, fill = type)) + 
    geom_histogram(position = 'dodge') 
  
  # TEST: 2-D densities 
  # 1. estimate 2-d prob density with MASS
  dens1 <- MASS::kde2d(adult_tidy$global_dist[adult_tidy$type == "pathway_dist"], adult_tidy$distance[adult_tidy$type == "pathway_dist"], n = 100)
  dens2 <- MASS::kde2d(adult_tidy$global_dist[adult_tidy$type == "random_dist"], adult_tidy$distance[adult_tidy$type == "random_dist"], n = 100)
  
  
  # 2. plot individual distributions 
  # create a tidy data frame from the density matrices
  dens_df <- expand.grid(x = dens1$x, y = dens1$y)
  dens_df$z <- as.vector(dens1$z)
  
  # create heatmap
  ggplot(dens_df, aes(x = x, y = y, fill = z)) +
    geom_tile() +    
    theme_minimal()
  
  
  # 3. Plot the difference 
  
  # Compute the difference in the 2D densities between the two distributions
  diff_dens <- dens1$z - dens2$z
  
  # create a tidy data frame from the density matrices
  dens_df <- expand.grid(x = dens1$x, y = dens1$y) # same axis for both
  dens_df$z <- as.vector(diff_dens)
  
  # create heatmap
  ggplot(dens_df, aes(x = x, y = y, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(dens_df$z)) +
    theme_minimal()
  
  ggsave(paste("./scripts/figures/", i, "_density_plt.pdf", sep=""))
  
}


