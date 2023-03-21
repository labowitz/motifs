library(Seurat)
library(stringr)

# Matrix of size 27998 x 92649

data_dir <- "data/raw_data/sc_data/gastrulation/"
sample_files <- list.files(path=data_dir, 
                           pattern = "GSE")
meta_files <- list.files(path=data_dir, 
                           pattern = "metadata")


# format is in 10x readable folders 
matrix_list = list() 
for( i in 1:length(sample_files))
  matrix_list[[i]] = Read10X(paste(data_dir, 
                                   sample_files[i], 
                                   sep=""))    

# Split a multi-page Excel file with different timepoints. 
sample_id = str_split(sample_files, 
                      '_', 
                      simplify=T)[,3]

# Read all csvs and save correct cell_id
# Cell column names will carry BC + sample id + stage
meta_list = list() 
for(i in 1:length(meta_files)){
  meta_list[[i]] = read.csv(paste(data_dir,
                                  meta_files[i], 
                                  sep=''),
                            
                            header = T)
  
  meta_list[[i]]$sample = sample_id[i]
  
  # define cell_id as BC + sample_id
  meta_list[[i]] <- meta_list[[i]] %>% 
    mutate(cell_id = paste(BC, 
                           sample, 
                           sep="_"))
}

# Merge all metadata files
meta_df = do.call(rbind, 
                  meta_list)

# Rename cells in matrix with BC + sample id + stage format
sample_id = str_split(sample_files,
                      '_', 
                      simplify  = T)[,3]

for(i in 1:length(matrix_list))
  colnames(matrix_list[[i]]) <- paste(colnames(matrix_list[[i]]), 
                                      sample_id[i], 
                                      sep="_")

# Merge all matrices
# All have the same number of genes so we can directly merge 
master_matrix = do.call(cbind, 
                        matrix_list)

# Make filter for Seurat 
which_cells = colnames(master_matrix)[which(colnames(master_matrix) %in% 
                                              meta_df$cell_id)]
master_matrix = master_matrix[,which_cells]

meta_df %>% 
  filter(cell_id %in% which_cells) -> meta_df

row.names(meta_df)<- meta_df$cell_id

# Seurat 
row.names(meta_df) <- meta_df$cell_id
tiss.embryo = CreateSeuratObject(counts = master_matrix, 
                                 min.cells = 2, 
                                 min.features = 800,
                                 meta.data = meta_df)

tiss.embryo = PercentageFeatureSet(tiss.embryo, 
                                   pattern = '^MT-', 
                                   col.name = 'percent.mt')

tiss.embryo <- tiss.embryo %>% 
  SCTransform() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 1.2)