library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)

tabula.path ="./data/raw_data/sc_data/tms_FACS/"

FACS_files = list.files(paste(tabula.path,"FACS" ,sep=""), full.names = TRUE)

raw.data.list = list()
for (file in FACS_files){
  #read each file (csv count matrices sorted by tissue)
  raw.data <- read.csv(file, row.names = 1)
  raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE)
  raw.data.list <- append(raw.data.list, raw.data)
}

raw.data <- do.call(cbind, raw.data.list)
cell_order_FACS <- order(colnames(raw.data))
raw.data = raw.data[,cell_order_FACS]

meta.data <- read.csv(paste(tabula.path,"/metadata_FACS.csv",sep=""))

plates <- str_split(colnames(raw.data),"[.]", simplify = TRUE)[,2]

rownames(meta.data) <- meta.data$plate.barcode

cell.meta.data <- meta.data[plates,]
rownames(cell.meta.data) <- colnames(raw.data) #lets rename the rows with the cell.ID

erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)

percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
sum.ercc<-Matrix::colSums(raw.data[erccs, ])
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,] #remove the ERCC sequences

x11();hist(percent.ercc)

total.reads = Matrix::colSums(raw.data)

tiss <- CreateSeuratObject(counts = raw.data, assay = "RNA")
tiss <- AddMetaData(object = tiss, cell.meta.data)
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
tiss<- AddMetaData(object = tiss, total.reads,col.name = 'total.reads')

ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(tiss), value = TRUE)
percent.ribo <- Matrix::colSums(tiss[ribo.genes, ])/Matrix::colSums(tiss)
tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")

#all of these percentages are for EACH cell
percent.Rn45s <- tiss@raw.data['Rn45s', ]/Matrix::colSums(tiss@raw.data)
tiss <- AddMetaData(object = tiss, metadata = percent.Rn45s, col.name = "percent.Rn45s")

FeatureScatter(object = tiss, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#filter cells that have less that 50000 reads, or that have more than 500 genes with zero reads
tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nReads"), low.thresholds = c(2, 50))

tiss <- NormalizeData(object = tiss, scale.factor = 1e6) #default normalization by Seurat
manual.norm = log(t( t(tiss@raw.data)/total.reads  ) * 1000000 +1)
tiss@data<-manual.norm

tiss <- ScaleData(object = tiss)
#this does not work for low number of low-expressed genes: see the PCA part
tiss <- FindVariableFeatures(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5,num.bin =4)
tiss <- FindVariableFeatures(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.2,num.bin =3,x.low.cutoff = 0.2,mean.function=ExpMean)


