source("./scripts/analysis/imports.R")

# Specify pathway name (as saved in pathway_df) and get genes.
pathway_name =  'Notch receptors, Dll ligands, and Fringe proteins'
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj = master_seurat)

# Thresholds for processing data.
min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.2 # tgfb minimum expression threshold for gene to be on
diverse_quantile = 0.9

subset = c("Mfng", "Rfng", "Lfng")

# Optimal number of pathway components (computed previously with silh. score)
optimal_k_pathway = 30

pipe_run <- quickPipeline(seurat_obj = master_seurat,
                          pathway_genes = pathway_genes,
                          k_final = optimal_k_pathway, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold
)

df <- gather(pipe_run$data_frame[,subset], gene, value, subset)
df$gene <- as.factor(df$gene)

ggplot(df, aes(x = value, group = gene, linetype = factor(gene))) + 
  stat_ecdf() + 
  theme_pubr() + 
  xlab("Normalized Expression") + 
  ylab("Empirical Cumulative Distribution") + 
  labs(linetype = "Fringe Gene") + 
  scale_linetype_manual(values = c(Lfng = "solid", Mfng = "dotted", Rfng = "dotdash"))

ggsave(paste(fig_dir, "Fringes_ECDF.pdf", sep=""), width = 7, height = 7)

