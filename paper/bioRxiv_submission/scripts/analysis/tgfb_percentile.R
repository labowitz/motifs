## Change these variables
filename = "tgfb" # how name is appended to files for pathway
pathway_index = 2 # Tgfb pathway index

## Directories for the files
source("./scripts/analysis/imports_percentile.R")

output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/percentile/", filename, sep = "")

pathway_name =  param_list[[pathway_index]][[1]] # Tgfb pathway name
min_genes_pathway = param_list[[pathway_index]][[2]] # Tgfb min. number of genes expressed
min_expr_threshold = 0.2 # Tgfb minimum expression threshold for gene to be on

optimal_k_pathway = 30 # optimal pathway components, computed from z-score

## Make a plot of k_opt profiles for cell types with pathway "ON
pipe_run <- quickPipeline(master_seurat = master_seurat, 
                          k_final = optimal_k_pathway, 
                          silent_plot = F, 
                          which_pathway = pathway_name, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold, 
                          which_profiles = 'both', 
                          save_pdf = T, 
                          pdf_file = paste(output_dir, "_profiles.pdf", sep = ""))

## Computing the optimal silhouette score
#silh_plt = silhouettePlot(which_pathway = pathway_name, 
#                          min_ON = min_genes_pathway, 
#                          min_expr = min_expr_threshold, 
#                          n_bootstraps = 100)

#g <- silhouette_zscore(silh_result = silh_plt,
#                       pathway_name = pathway_name,
#                       x_offset = 6,
#                       max.y = 0.43,
#                       min.y =0.17 # adjust axis parameters 
#)
#g
#ggsave(paste(output_dir, "_silh_z-score.pdf", sep = ""), width = 7, height = 7)

## Running fullControlPathway
control_res  = fullControlPathway(this_pathway = pathway_name,
                                  k_pathway = optimal_k_pathway , 
                                  filter_pathway = 'both', # adult + devel 
                                  this_seurat = master_seurat, # seurat object
                                  null_list = hvg_genes, #list of highly variable genes 
                                  n_samples = 100, filter_manual = T,
                                  min_expr_gene = min_expr_threshold, 
                                  min_genes_ON = min_genes_pathway, # default > 2 genes ON 
                                  n_pcs_global = 100, # how many PCs to use
                                  embedding_matrix = pca_proj # PCA embedding for diversity 
)


## Plotting
# Violin plots
violin_plot = violin_diversity(control_res)
violin_plot
ggsave(paste(output_dir, "_violin_plot.pdf", sep = ""), width = 7, height = 7)

# ECDF plots
ecdf_plot = ecdf_diversity(control_res)
ecdf_plot
ggsave(paste(output_dir, "_ecdf_plot.pdf", sep = ""), width = 7, height = 7)

# Diversity rank plots
rank_plot <- rank_diversity(
  which_pathway = pathway_name, 
  k_motifs = optimal_k_pathway, 
  min_expression = min_expr_threshold, 
  min_genes_pathway = min_genes_pathway,
  embedding_matrix = pca_proj
)
rank_plot
ggsave(paste(output_dir, "_rank_plot.pdf", sep = ""), width = 7, height = 7)

# Plot global dendrogram with pathway states
glob_dendr <- globalDendrogramPlot2(
  control_res = control_res, 
  seurat_obj = master_seurat, 
  use_pca = T, 
  npcs = 1:20,
  clust_method = 'ward.D2',
  dist_metric ='cosine',
  save_dir = output_dir
)

## Plot the most diverse profiles, aka the motifs
# 1. New: Select all profiles 
diverse_df <- control_res$profiles
# 2. tidy data.frame to average gene expression within a pathway profile 
diverse_df %>% pivot_longer(cols = genesPathway(pathway_name), 
                            names_to = 'gene', 
                            values_to = 'expression') %>% 
  select(cell_id, cell_ontology_class, Tissue, Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy 

# 3. average profile 
diverse_tidy %>% group_by(class_label,gene) %>% summarise(mean_expr = mean(expression), 
                                                          rank = mean(rank),diversity = mean(diversity), 
                                                          cell_types = mean(n)) -> diverse_tidy

# 4. wide format 
diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types), names_from=gene,values_from=mean_expr) %>% tibble::column_to_rownames('class_label')-> diverse_mat

# 5. We do the filtering here either for motifs or for non-diverse profiles 
control_res$diversity %>% 
  dplyr::filter(type=='transcriptome') %>% # choose the null model 
  pull(d) %>% quantile(0.90) -> divers_quantile

diverse_mat %>% dplyr::filter(diversity>divers_quantile) -> diverse_mat 

pdf(paste(output_dir, "_motif_profiles.pdf", sep = ""))
motif_heatmap <- superheat(diverse_mat[,genesPathway(pathway_name)],
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
dev.off()

# Motif profile cell type distribution
diverse_df %>% dplyr::filter(class_label %in% row.names(diverse_mat) ) %>% 
  select(Tissue, class_label ) %>% dplyr::filter(!Tissue %in% c('mat','endoderm')) %>% group_by(class_label,Tissue)%>% 
  count  %>% 
  pivot_wider(id_cols=class_label,names_from = Tissue,
              values_from= n,values_fill = 0) -> tissue_distribution


tissue_pal<-colorRampPalette(brewer.pal(n = 9, name = 'PuRd'))
x = tissue_distribution %>% ungroup %>% select(-class_label) %>% as.matrix() 
row.names(x) <- tissue_distribution$class_label
# make tissue names pretty 
colnames(x)<-str_replace(string = colnames(x),pattern ="_",replacement = ' ') %>% firstup() 

x <- x[,sort(colnames(x))]

pheatmap(sqrt(x), treeheight_row = 20,treeheight_col = 20,
         clustering_method = 'ward.D2',col = tissue_pal(100),
         cluster_rows = F,
         cluster_cols = F,
         fontsize =12,angle_col = 45,
         filename = paste("./scripts/analysis/outputs/tgfb_analysis/", filename, "_tissue_distribution.pdf", sep = ""),
         height =4, width = 6) # keep this size to make it similar to the motif heatmap 

# Private profiles
# 1. New: Select all profiles 
diverse_df <- control_res$profiles
# 2. tidy data.frame to average gene expression within a pathway profile 
diverse_df %>% pivot_longer(cols = genesPathway(pathway_name), 
                            names_to = 'gene', 
                            values_to = 'expression') %>% 
  select(cell_id, cell_ontology_class, Tissue, Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy 

# 3. average profile 
diverse_tidy %>% group_by(class_label,gene) %>% summarise(mean_expr = mean(expression), 
                                                          rank = mean(rank),diversity = mean(diversity), 
                                                          cell_types = mean(n)) -> diverse_tidy

# 4. wide format 
diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types), names_from=gene,values_from=mean_expr) %>% tibble::column_to_rownames('class_label')-> diverse_mat

# 5. We do the filtering here either for motifs or for non-diverse profiles 
control_res$diversity %>% 
  dplyr::filter(type=='transcriptome') %>% # choose the null model 
  pull(d) %>% quantile(0.40) -> divers_quantile

diverse_mat %>% dplyr::filter(diversity<divers_quantile) -> diverse_mat 

pdf(paste(output_dir, "_private_profiles.pdf", sep = ""))
motif_heatmap <- superheat(diverse_mat[,genesPathway(pathway_name)],
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
dev.off()


# Save the data.frame
glob_dendr[[2]] %>% group_by(class_label, cell_ontology_class) %>% summarise(age = unique(age)) %>% write.csv(paste(output_dir, "_profile_cell_ontology_class.csv", sep = ""))