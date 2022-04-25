library(forcats)
source("./scripts/analysis/imports.R")
output_dir <- "./scripts/figures/"

# Prepare df
df <- master_seurat@meta.data %>% dplyr::select( cell_ontology_class,age,dataset,Cell_class,Tissue ) %>% mutate(full_celltype = paste(cell_ontology_class, age, Tissue ) ) %>% ungroup() %>% select( dataset, full_celltype, Cell_class) %>% group_by(dataset, Cell_class)  %>% unique() %>% count()

# Reorder df in order of timing
dataset_order <- rev(c("E5.5_Nowotschin", "E6.5_8.5_Chan", "E6.5_8.5_Marioni", "E9.5_11.5_Tang", "Forelimb_E10.5_15.0", "1m 10x", "3m 10x", "FACS 3m", "18m 10x", "FACS 18m", "21m 10x", "24m 10x", "FACS 24m", "30m 10x"))
df <- df %>% arrange(match(dataset, dataset_order), Cell_class, n)

# Make plot
p <- ggplot(df, aes(x=n, y=fct_inorder(dataset), fill=Cell_class)) + 
  geom_bar(stat="identity") + scale_fill_manual(values=colors_1206$Cell_class) + 
  theme_bw()
ggsave(paste0(output_dir, "Figure_2A.pdf", sep=""))
