library(Seurat)
library(dplyr)
library(stringr)

wnt_ligands = c('Wnt1', 'Wnt2', 'Wnt2b', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5a',
        'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Wnt8b', 'Wnt9a', 'Wnt9b',
        'Wnt10a', 'Wnt10b', 'Wnt11', 'Wnt16')

wnt_receptors = c('Lrp5', 'Lrp6', 'Fzd1', 'Fzd2', 'Fzd3',
        'Fzd4', 'Fzd5', 'Fzd6', 'Fzd7', 'Fzd8', 'Fzd9', 'Fzd10')

bmp_ligands = c('Bmp2', 'Bmp3', 'Bmp4', 'Bmp5', 'Bmp6', 'Bmp7', 'Bmp8a',
        'Bmp10', 'Bmp11', 'Bmp15', 'Gdf6', 'Gdf7', 'Gdf5', 'Gdf10',
        'Gdf11')

bmp_receptors = c("Bmpr1a","Bmpr1b", "Acvr1", "Acvrl1", "Acvr1b", "Tgfbr1","Acvr1c",
        "Acvr2a", "Acvr2b", "Bmpr2", "Tgfbr2")

notch = c("Dll1", "Dll3","Dll4", "Jag1", "Jag2", "Notch1", "Notch2", 
         "Notch3", "Notch4", "Mfng", "Rfng", "Lfng")

eph_receptors = c('Epha1', 'Epha2', 'Epha3', 'Epha4', 'Epha5', 'Epha7', 'Ephb1', 
        'Ephb2', 'Ephb3', 'Ephb4', 'Ephb6')

eph_ligands = c('Efna1', 'Efna2', 'Efna3', 'Efna4', 'Efna5', 'Efnb1', 
        'Efnb2', 'Efnb3')

fgfr = c('Fgf1', 'Fgf10', 'Fgf11', 'Fgf12', 'Fgf13', 'Fgf14', 'Fgf16', 'Fgf18', 
        'Fgf2', 'Fgf7', 'Fgf9', 'Fgfbp1', 'Fgfbp3', 'Fgfr1', 'Fgfr1op', 'Fgfr1op2', 
        'Fgfr2', 'Fgfr3', 'Fgfr4', 'Fgfrl1')

splice_srsf = c('Srsf1', 'Srsf10', 'Srsf11', 'Srsf12', 'Srsf2', 'Srsf3', 'Srsf4', 
               'Srsf5', 'Srsf6', 'Srsf7', 'Srsf9')

lpa = c('Lpar1', 'Lpar2', 'Lpar3', 'Lpar4', 'Lpar5', 'Lpar6')

all_genes = c(wnt_ligands, wnt_receptors, bmp_ligands, bmp_receptors, notch,
  eph_receptors, eph_ligands, fgfr, splice_srsf, lpa)

all_genes <- toupper(all_genes)

all_pathways <- readRDS("data/processed_data/all_pathways.RDS")

data_dir <- "data/raw_data/sc_data/gut_endoderm/"
sample_files <- list.files(path=data_dir, 
                           pattern = "_counts.csv")
meta_data_file <- list.files(path=data_dir, 
                             pattern = "_metadata.csv")

meta_author <- read.csv(paste(data_dir, 
                            meta_data_file, 
                            sep=""))

cell_ids = c('E5.5_1', 'E5.5_2', 'E5.5_3')
sample_id = stringr::str_split(
  stringr::str_split(sample_files, 
                      "counts", 
                      simplify = T)[,1],
                      "_", n=2, 
  simplify = T)[,2]

# Assuming all files respect the same format 
matrix_list = list() 
for( i in 1:length(sample_files)){
  
  a = read.csv(paste(data_dir, 
                     sample_files[i], 
                     sep=""), 
               header = T)
  
  a[,-1] %>% as.matrix() %>% t() -> x_mat
  
  colnames(x_mat) = paste(a$X, 
                          cell_ids[i], 
                          sep = "_")
  
  barcodes = str_split(colnames(x_mat), 
                       "\\'", 
                       simplify = T)[,2]
  
  cell_id = paste(sample_id[i],
                  barcodes, 
                  sep ='')
  
  colnames(x_mat) <- cell_id
  
  matrix_list[[i]] = x_mat
  
  # Fix column names to match author's meta data 
}


imputeAbsetGenes <- function(x_mat, 
                             which_genes){
  
  #which_genes = all_pathways$gene[all_pathways$pathway=='BMP 7receptors']
  
  #x_mat = matrix_list[[1]]
  not_present = which_genes[which(!toupper(which_genes) %in% row.names(x_mat))    ]
  
  fill_in = matrix(0, length(not_present), dim(x_mat)[2])
  row.names(fill_in)<-not_present
  colnames(fill_in) <- colnames(x_mat)
  
  return(rbind(x_mat, fill_in))
}


matrix_list_fil = lapply(matrix_list, 
                         imputeAbsetGenes, 
                         all_pathways$gene %>% 
                           unique() %>% 
                           toupper())
genes_list = lapply(matrix_list_fil, row.names)
shared_genes = Reduce(intersect, genes_list)

master_matrix = matrix_list_fil[[1]][shared_genes,]
for( i in 2:length(matrix_list_fil)){
  
  master_matrix = cbind(master_matrix, 
                        matrix_list_fil[[i]][shared_genes,])
}

meta.data = data.frame(cell = colnames(master_matrix), 
                       state = str_split(colnames(master_matrix), 
                                         '_', 
                                         simplify = T)[,2]) 
row.names(meta.data) <- meta.data$cell

meta.data$batch = paste(meta.data$state,
                        str_split(meta.data$cell, 
                                  '_', 
                                  simplify =T)[,3], 
                        sep ="_")

# Create seurat objects 
tiss.embryo = list() 
for(i in 1:length(matrix_list)){
  which_cells = colnames(matrix_list[[i]])[which(colnames(matrix_list[[i]]) %in% meta_author$index)]
  sub_meta = meta_author %>% 
    filter(index %in% which_cells)    
  row.names(sub_meta) <- sub_meta$index
  tiss.embryo[[i]] = CreateSeuratObject(counts= matrix_list[[i]][,which_cells],
                                        meta.data = sub_meta)
}

# Apply SCT transform to all samples in the list 
for(i in 1:length(tiss.embryo)){
  tiss.embryo[[i]] <- SCTransform(tiss.embryo[[i]], 
                                  verbose = F)
}

# Change the memory parameters for parallel computing (future)
# this dataset requires 5.6Gb so I set it to 6Gb
options(future.globals.maxSize = 6000 * 1024^2)

# select features for downstream analysis 

embryo.features <- SelectIntegrationFeatures(object.list = tiss.embryo, 
                                             nfeatures = 3000)
tiss.embryo <- PrepSCTIntegration(object.list = tiss.embryo, 
                                  anchor.features = embryo.features, 
                                  verbose = FALSE)

# Manually include all signaling genes if not present 
embryo.features = c(embryo.features, 
                    all_genes[which(!all_genes %in% embryo.features)])
# Integrate

embryo.anchors <- FindIntegrationAnchors(object.list = tiss.embryo, 
                                         normalization.method = "SCT", 
                                         anchor.features = embryo.features, 
                                         verbose = FALSE)

embryo.integrated <- IntegrateData(anchorset = embryo.anchors, 
                                   normalization.method = "SCT", 
                                   verbose = FALSE)


# From here, normal pipeline
embryo.integrated <- RunPCA(embryo.integrated, 
                            verbose = FALSE)
embryo.integrated <- RunUMAP(embryo.integrated, 
                             dims = 1:30)


embryo.integrated %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters( resolution = 1.5) -> embryo.integrated

