source("./scripts/analysis/imports.R")

adult_subset <- master_seurat[,master_seurat$dataset %>% str_detect("10x|FACS")]

n_pcs = 30

global_coords = Embeddings(master_seurat, reduction='pca')
global_coords = global_coords[adult_subset$cell_id, 1:n_pcs] # already set row names as cell id
rownames(global_coords) <- adult_subset$cell_ontology_class

cell_type_dist_df <- data.frame()

for (i in unique(adult_subset$cell_ontology_class)){
  
  if ((adult_subset$cell_ontology_class == i) %>% sum() > 1){
    
    global_dist = global_coords[rownames(global_coords) == i,]
    
    centroids = global_dist %>% colMeans()
    
    global_dist = (global_dist - centroids)^2 %>% rowSums() %>% sqrt()
    
    global_dist_df <- data.frame(i, 
                                 as.numeric(global_dist))
    
    colnames(global_dist_df) <- c("cell_type", "global_dist")
    
    cell_type_dist_df <- rbind(cell_type_dist_df, global_dist_df)
    
  }
}

p <- ggplot(cell_type_dist_df, aes(cell_type, global_dist)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
