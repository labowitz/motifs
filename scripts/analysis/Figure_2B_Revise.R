source("./scripts/analysis/imports.R")

## Figure 2B bottom

load("./data/processed_data/tiss_SCT_Jan2020_processed.rdata")

# Add cell class annotations
df <- data.frame(seurat_clusters = 0:81, cell_ontology_class = rep("0", 82))

for(i in 1:82) {
  df$cell_ontology_class[[i]] = tiss.norm@meta.data %>% filter(seurat_clusters == (i-1)) %$% cell_ontology_class %>% table %>% which.max %>% names
}

imp_df <- read.csv("./data/processed_data/countMat_FACS_Tabula3m_meta.csv")

df$Cell_class <- imp_df$cell_class[match(df$seurat_clusters, df$seurat_clusters)]

identical(df$Cell_class, imp_df$cell_class)

tiss.norm@meta.data$Cell_class <- df$Cell_class[match(tiss.norm@meta.data$seurat_clusters, df$seurat_clusters)]

# Plot the UMAP
DimPlot(tiss.norm, 
        reduction = "umap",
        group.by = "Cell_class", # some column in the meta.data of the object
        cols = colors_1206$Cell_class
)

ggsave(paste0(output_dir, "Figure_2B_left_bottom.pdf", sep=""))