source("./scripts/analysis/imports.R")
fig_dir = "./scripts/figures/"

## Integrated atlas
DimPlot(master_seurat, 
        reduction = "umap",
        group.by = "Cell_class", # some column in the meta.data of the object
        cols = colors_1206$Cell_class
)

ggsave(paste0(fig_dir, "Figure_2B_main.pdf", sep=""))