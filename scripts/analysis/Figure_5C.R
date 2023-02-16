source("./scripts/analysis/imports_new.R")
source("./scripts/analysis/Functions_Figure_5.R") # Script with Fig. 5 functions

#silh_res_dir = "data/processed_data/Silhouette_PathBank/"
#dispersion_dir = "data/processed_data/Dispersion_PathBank/"

silh_res_dir = "scripts/figures/peak_analysis/silhouette_res/silh_rds/"
dispersion_dir = "scripts/figures/peak_analysis/dispersion/"

# Parameters to show
min_expr_threshold <- 0.3
min_genes_pathway <- 2

## Saved silhouette scores
silh_files <- list.files(path = silh_res_dir,
                               pattern = ".RDS")

df_kvals <- peak_width_scores(use_percentile = 0.95, 
                              smooth_window = 3, 
                              pct_opt = 0.9, 
                              peak_threshold = 0.1,
                              silh_res_dir = silh_res_dir,
                              silh_files = silh_files,
                              min_expr = min_expr)

## Running code
df_kvals$width %>% hist(breaks = 20 , 
                        col = 'lightgray', 
                        border = 'black', 
                        xlab = 'Width',
                        main = "Peak sharpness", 
                        prob = T)

computeDispersions(df_kvals = df_kvals,
                   overwrite_files = T,
                   use_min_k = F,
                   silh_res_dir = silh_res_dir,
                   dispersion_dir = dispersion_dir)

dispersion_stats <- parseDispersion(pathway_list_dispersion = pathway_list_dispersion[2:length(pathway_list_dispersion)],
                                    dispersion_dir = dispersion_dir)

# 3. Master data.frame 
res = plot_dispersion_distribution(dispersion_stats ,
                                   labels_idx = c(9,10,12,13,14,1,3), 
                                   min_silh = 0.3,
                                   min_n_cluster = 2,  
                                   df_kvals = df_kvals, 
                                   use_mean = F,
                                   include_metabolic = F, 
                                   use_palette = 'Paired')

res[[1]] + 
  theme(text = element_text(size = 20)) +
  ylab('N clusters / n genes') + 
  xlab('Dispersion')

master_df <- res[[2]] # will now include all the summary statistics for all pathways
master_df <- master_df %>% mutate(short_label = str_trunc(pathway_name, 10))

fig5_scatter(master_df = master_df, max_width = 0.35)
