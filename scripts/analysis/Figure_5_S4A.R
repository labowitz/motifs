source("./scripts/analysis/imports.R")
source("./scripts/analysis/Figure_5_Functions.R") # Script with Fig. 5 functions

# Ran the silhoeutte score and dispersion computations earlier and stored here
silh_res_dir = "./scripts/figures/peak_analysis/silhouette_res/silh_rds/alejo_res/"
dispersion_dir = "./scripts/figures/peak_analysis/dispersion/"

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
min_genes_pathway = 2
min_expr_threshold = 0.3
diverse_quantile = 0.9

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
pathway_list_dispersion <- full_pathway_list[saved_idx]

df_kvals <- peak_width_scores(use_percentile = 0.95, 
                              smooth_window = 3, 
                              pct_opt = 0.9, 
                              peak_threshold = 0.1,
                              silh_res_dir = silh_res_dir,
                              silh_files = saved_files,
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
                   dispersion_dir = dispersion_dir,
                   silh_files = saved_files)

dispersion_stats <- parseDispersion(pathway_list_dispersion = pathway_list_dispersion[2:length(pathway_list_dispersion)],
                                    dispersion_dir = dispersion_dir)

# 3. Master data.frame 
res = plot_dispersion_distribution(dispersion_stats,
                                   labels_idx = c(9,10,12,13,14,1,3), 
                                   min_silh = 0.3,
                                   min_n_cluster = 2,  
                                   df_kvals = df_kvals, 
                                   use_mean = F,
                                   include_metabolic = F, 
                                   use_palette = 'Paired',
                                   silh_files = saved_files)

res[[1]] + 
  theme(text = element_text(size = 20)) +
  ylab('N clusters / n genes') + 
  xlab('Dispersion')

master_df <- res[[2]] # will now include all the summary statistics for all pathways
master_df <- master_df %>% mutate(short_label = str_trunc(pathway_name, 10))

fig5_scatter(master_df = master_df, max_width = 0.35)
