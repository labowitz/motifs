source("./scripts/analysis/imports_new.R")

## Change these variables
filename = "tgfb_trial" # how name is appended to files for pathway

## Directories for the files
output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")
fig_dir = "./scripts/figures/tgfb_trial/"

pathway_name =  "Bmp_Tgfb" # tgfb pathway name

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
optimal_k_pathway = 30 # optimal pathway components, computed from z-score

## Generate the dataframe with the normalized gene counts of profiles with pathway ON
sat_val = 0.99
diverse_quantile = 0.9
seurat_obj = master_seurat
min_expr = min_expr_threshold
min_genes_on = min_genes_pathway
k_final = optimal_k_pathway

res_list = quickPipeline(pathway_genes = pathway_genes,
                         seurat_obj = seurat_obj,
                         k_final = k_final,
                         min_genes_on=min_genes_on,
                         min_expr=min_expr
)

## Generate the correlation matrix.
M_corr = as.data.frame(cor(res_list$data_frame %>% dplyr::select(pathway_genes)))
colnames(M_corr) = pathway_genes
rownames(M_corr) = pathway_genes

## Correlation plot for all profiles.
pdf(file = paste0(fig_dir, "bmpr_corrplot_all.pdf", sep=""), useDingbats=FALSE)
corrplot(as.matrix(M_corr), type = "upper")
dev.off()

## Correlation plot for motif profiles
# Get the motifs
div_res <- diverseFilt(control_res = control_res,
                       pathway_genes = pathway_genes,
                       diverse_quantile = diverse_quantile,
                       type="motif"
)

M_corr = as.data.frame(cor(res_list$matrix[rownames(div_res$diverse_mat),] %>% dplyr::select(pathway_genes)))
colnames(M_corr) = pathway_genes
rownames(M_corr) = pathway_genes

pdf(file = paste0(fig_dir, "bmpr_corrplot_motifs.pdf", sep=""), useDingbats=FALSE)
corrplot(as.matrix(M_corr), type = "upper")
dev.off()

## Correlation plot for private profiles
# Get the motifs
div_res <- diverseFilt(control_res = control_res,
                       pathway_genes = pathway_genes,
                       diverse_quantile = diverse_quantile,
                       type="private"
)

M_corr = as.data.frame(cor(res_list$matrix[rownames(div_res$diverse_mat),] %>% dplyr::select(pathway_genes)))
colnames(M_corr) = pathway_genes
rownames(M_corr) = pathway_genes

pdf(file = paste0(fig_dir, "bmpr_corrplot_private.pdf", sep=""), useDingbats=FALSE)
corrplot(as.matrix(M_corr), type = "upper")
dev.off()