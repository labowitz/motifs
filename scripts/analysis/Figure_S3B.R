library(corrplot)
source("./scripts/analysis/imports_new.R")

fig_dir = "./scripts/figures/"

pathway_name <- "Tgf-beta family receptors" # tgfb pathway name
pathway_genes <- genesPathway(pathway_name = pathway_name,
                              pathway_df = pathway_df)

## Generate the dataframe with the normalized gene counts of profiles with pathway ON
sat_val = 0.99
diverse_quantile = 0.9
seurat_obj = master_seurat
min_expr_threshold = 0.2
min_genes_pathway = 2
optimal_k_pathway = 30

res_list = quickPipeline(pathway_genes = pathway_genes,
                         seurat_obj = seurat_obj,
                         k_final = optimal_k_pathway,
                         min_genes_on=min_genes_pathway,
                         min_expr=min_expr_threshold
)

control_res  = fullControlPathway(pathway_genes = pathway_genes,
                                  k_final = optimal_k_pathway,
                                  seurat_obj = master_seurat, # seurat object
                                  n_samples = 100, 
                                  filter_manual = T,
                                  min_genes_on = min_genes_pathway, 
                                  min_expr = min_expr_threshold, 
                                  n_pcs = 100, # how many PCs to use
                                  manual_embedding = pca_proj, # PCA embedding for diversity 
                                  dist_metric = "euclidean"
)

## Generate the correlation matrix.
M_corr = as.data.frame(cor(res_list$data_frame %>% 
                             dplyr::select(pathway_genes)))
colnames(M_corr) = pathway_genes
rownames(M_corr) = pathway_genes

## Correlation plot for all profiles.
pdf(file = paste0(fig_dir, 
                  "Figure_S3B_all.pdf", 
                  sep=""), 
    useDingbats=FALSE)
corrplot(as.matrix(M_corr), type = "upper")
dev.off()

## Correlation plot for motif profiles
# Get the motifs
div_res <- diverseFilt(control_res = control_res,
                       pathway_genes = pathway_genes,
                       diverse_quantile = diverse_quantile,
                       type="motif"
)

M_corr = as.data.frame(cor(res_list$matrix[rownames(div_res$diverse_mat),] %>% 
                             dplyr::select(pathway_genes)))
colnames(M_corr) = pathway_genes
rownames(M_corr) = pathway_genes

pdf(file = paste0(fig_dir, 
                  "Figure_S3B_motif.pdf", sep=""), 
    useDingbats=FALSE)
corrplot(as.matrix(M_corr), type = "upper")
dev.off()

## Correlation plot for private profiles
# Get the motifs
div_res <- diverseFilt(control_res = control_res,
                       pathway_genes = pathway_genes,
                       diverse_quantile = diverse_quantile,
                       type="private"
)

M_corr = as.data.frame(cor(res_list$matrix[rownames(div_res$diverse_mat),] %>% 
                             dplyr::select(pathway_genes)))
colnames(M_corr) = pathway_genes
rownames(M_corr) = pathway_genes

pdf(file = paste0(fig_dir, 
                  "Figure_S3B_private.pdf", 
                  sep=""), 
    useDingbats=FALSE)
corrplot(as.matrix(M_corr), type = "upper")
dev.off()