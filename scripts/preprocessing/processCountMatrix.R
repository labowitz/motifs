
# ---
# title: "Tabula Muris: All FACS analysis"
# output: html_notebook
# ---

# Preprocessing
#
# Load the requisite packages and some additional helper functions.
#
# ```{r}
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(readr)
library(here)


#get the annotated genes (manually curated lists)
pathway.genes<-function(pathway ="bmp"){
  bmp.receptors<-c("Bmpr1a","Bmpr1b","Acvr1","Acvrl1","Acvr1b","Tgfbr1","Acvr1c","Acvr2a","Acvr2b","Bmpr2","Tgfbr2")
  bmp.ligands<-c("Bmp2","Bmp3","Bmp4","Bmp5","Bmp6","Bmp7",
              "Bmp8a","Gdf3","Gdf9","Gdf10","Gdf11","Gdf15")
  bmp.smads<-c("Smad1" ,"Smad2" ,"Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9")

  notch.all<-c(
  "Dll1",
  "Dll3",
  "Dll4",
  "Dtx1",
  "Jag1",
  "Jag2",
  "Adam10",
  "Psen1",
  "Psen2",
  "Psenen",
  "Notch1",
  "Notch2",
  "Notch3",
  "Notch4")


    if(pathway =="bmp"){
      genes.plot = c(bmp.receptors,bmp.ligands,bmp.smads)
    }else if(pathway=="notch"){
      genes.plot = notch.all
    }

    return (genes.plot)
}
#optional:

tabula.path ="./data/raw_data/tabula_muris_facs/"

FACS_files = list.files(paste(tabula.path,"/FACS" ,sep=""), full.names = TRUE)

raw.data.list = list()
for (file in FACS_files){
  #read each file (csv count matrices sorted by tissue)
  # FACS_files[1] ="/home/agranado/MEGA/Caltech/rnaseq/tabula-muris/00_data_ingest/00_facs_raw_data/FACS/Aorta-counts.csv"
  raw.data <- read.csv(file, row.names = 1)
  raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE)
  raw.data.list <- append(raw.data.list, raw.data)
}

#dim(raw.data.list[[2]])
#  [1] 23433  1638 # 1638 cells in this tissue, 23433 genes for ALL tissues]
# join the matrix with do.call cbind
raw.data <- do.call(cbind, raw.data.list)

  cell_order_FACS <- order(colnames(raw.data))
  raw.data = raw.data[,cell_order_FACS]

meta.data <- read.csv(paste(tabula.path,"metadata_FACS.csv",sep=""))
#we extract information from the colnames (ID of individual cells)

# meta.data includes the information for all plates

# > head(meta.data)
#           plate.barcode mouse.id  tissue subtissue FACS.selection mouse.sex
# D041914         D041914    3_8_M Bladder                 Multiple         M
# D042253         D042253    3_9_M Bladder                 Multiple         M
# MAA000487     MAA000487   3_10_M Bladder                 Multiple         M
# B000610         B000610   3_56_F Bladder                 Multiple         F
# B002764         B002764   3_38_F Bladder                 Multiple         F
# B002771         B002771   3_39_F Bladder                 Multiple         F

plates <- str_split(colnames(raw.data),"[.]", simplify = TRUE)[,2]

rownames(meta.data) <- meta.data$plate.barcode
#let's extract the plate information for all cells,
#plates is an array with the plate.id for each cell in the raw.data matrix
cell.meta.data <- meta.data[plates,]
rownames(cell.meta.data) <- colnames(raw.data) #lets rename the rows with the cell.ID
#because

# Find ERCC's, compute the percent ERCC, and drop them from the raw data.
# External RNA Controls Consortium
# The controls consist of a set of unlabeled, polyadenylated transcripts designed to be
# added to an RNA analysis experiment after sample isolation, in order to measure against defined performance criteria.

#look for more about ERCC here http://tools.thermofisher.com/content/sfs/manuals/cms_086340.pdf
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
# there are genes (rows) that come from the ERCC control, so we can find them and calibrate the quantifitation
# there should be 92 ERCC transcripts

#percent is the ratio, for each cell, between the sum of ERCC detection divided by the total count
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
sum.ercc<-Matrix::colSums(raw.data[erccs, ])
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,] #remove the ERCC sequences

x11();hist(percent.ercc)
 # # # FILTER
raw.data = raw.data[,percent.ercc<0.15] #Here is were the distribution seems to change
#filter on ERCC % which is a measure of RNA content. High ERCC means most of the sample comes from spike-in sequences

# # # # #HERE WE DO BMP-specific analysis
# lets filter manually before selecting bmp genes:
# FILTER CELLS
high.count.cells = Matrix::colSums(raw.data)>50000 #recommendation for quality control (Sisi)
raw.data  = raw.data[,high.count.cells]




######lets save the total count of reads,
# Seurat normalizes the data assuming that genes in the matrix are transcriptome
# therefore  gene_i/sum(gene)  is biased when taking only a subset of genes and we get weird things

  total.reads = Matrix::colSums(raw.data)


 # bmp.all = pathway.genes(pathway = "bmp")
 # bmp.indexes=match(bmp.all,rownames(raw.data))

 # raw.data <-raw.data[bmp.indexes,]


  #counts.cells = Matrix::colSums(raw.data)
 # raw.data = raw.data[,counts.cells>length(bmp.all)*2] #keep cells that have at least this N of counts


#  raw.logical = raw.data>0
#  genes.at.least.one = Matrix::colSums(raw.logical)
#  raw.data= raw.data[,genes.at.least.one>=3] #ten percent of genes in the pathway (31 genes)


  #explore the raw.data  :
  #would be good to filter based on CV but we can't do that yet...TO DO
#  cv.bmp <-apply(raw.data,1,sd)/apply(raw.data,1,mean)

#take from meta.data only cells that remain in the matrix
cell.meta.filter = cell.meta.data[colnames(raw.data), ]


#take the total read count only for the cells that remain
total.reads = total.reads[colnames(raw.data)]

#last step: filter duplicated cells (with exact same count profiles)
duplicated.cells = duplicated(t(as.matrix(raw.data)))
raw.data <- as.matrix(raw.data)
raw.data=raw.data[,!duplicated.cells]
#remove the duplicate from the metadata
cell.meta.filter=cell.meta.filter[colnames(raw.data),]
total.reads = total.reads[colnames(raw.data)]
#-------------------------------------------------------------------------------

# # # # ## SEURAT OBJECT
tiss <- CreateSeuratObject(raw.data = raw.data)

tiss <- AddMetaData(object = tiss, cell.meta.data)
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
tiss<- AddMetaData(object = tiss, total.reads,col.name = 'total.reads')
# Change default name for sums of counts from nUMI to nReads
colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads' # this is not UMI data so Seurat calculates only the number of reads

# > tiss
# An object of class Seurat
# 23341 features across 53760 samples within 1 assay
# Active assay: RNA (23341 features)

  ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
  percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
  tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")

#re write for Seurat v3:
#ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(tiss), value = TRUE)
#percent.ribo <- Matrix::colSums(tiss[ribo.genes, ])/Matrix::colSums(tiss)
#tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")


# #all of these percentages are for EACH cell
# percent.Rn45s <- tiss@raw.data['Rn45s', ]/Matrix::colSums(tiss@raw.data)
# tiss <- AddMetaData(object = tiss, metadata = percent.Rn45s, col.name = "percent.Rn45s")

# Analysis

#A sanity check: genes per cell vs reads per cell.


GenePlot(object = tiss, gene1 = "nReads", gene2 = "nGene", use.raw=T)

# # # # #Filter out cells with few reads and few genes.

#filter cells that have less that 50000 reads, or that have more than 500 genes with zero reads
tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nReads"), low.thresholds = c(2, 50))


# # #Normalize the data, then center and scale.
# Equivalent to
# # # r = t(t(raw.data)/Matrix::colSums(raw.data))
#     log(r[1:3,1:3] * 1000000 + 1)

#make Seurat think we normalized using their method...
tiss <- NormalizeData(object = tiss, scale.factor = 1e6) #default normalization by Seurat

#lets try to force normalization by hand :
manual.norm = log(t( t(tiss@raw.data)/total.reads  ) * 1000000 +1)
tiss@data<-manual.norm

tiss <- ScaleData(object = tiss)
#this does not work for low number of low-expressed genes: see the PCA part
tiss <- FindVariableFeatures(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5,num.bin =4)
tiss <- FindVariableFeatures(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.2,num.bin =3,x.low.cutoff = 0.2,mean.function=ExpMean)

#for BMP we see 11 variable genes

# # # #Run Principal Component Analysis.
#The RunPCA() function performs the PCA on genes in the @var.genes slot
#by default and this can be changed using the pc.genes parameter.
tiss <- RunPCA(object = tiss, do.print = FALSE, pcs.compute = 13)
#use all genes of the pathway, instead of only the "variable genes"
tiss <- RunPCA(object = tiss, pc.genes =bmp.genes, do.print = FALSE, pcs.compute = 25)


tiss <- ProjectPCA(object = tiss, do.print = FALSE)


#```{r, echo=FALSE, fig.height=4, fig.width=8}
PCHeatmap(object = tiss, pc.use = 1:3, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, num.genes = 8)
#```

# Later on (in FindClusters and TSNE) you will pick a number of principal components to use. This has the effect of keeping the major directions of variation in the data and, ideally, supressing noise. There is no correct answer to the number to use, but a decent rule of thumb is to go until the plot plateaus.


PCElbowPlot(object = tiss, num.pc = 20)

#
# Choose the number of principal components to use.
# ```{r}
# Set number of principal components.
n.pcs = 10
# ```
#
# The clustering is performed based on a nearest neighbors graph.
#Cells that have similar expression will be joined together.
#The Louvain algorithm looks for groups of cells with high modularity--more connections within the group than between groups.
#The resolution parameter determines the scale...higher resolution will give more clusters, lower resolution will give fewer.
#
# For the top-level clustering, aim to under-cluster instead of over-cluster. It will be easy to subset groups and further analyze them below.
#
# ```{r}
# Set resolution
#more resolution means more clusters:
#Seurat recommends between 0.6 and 1.2
#Tabula Muris papers uses 1 for 40k cells
res.used <- 1.2

tiss <- FindClusters(object = tiss, reduction.type = "pca", dims.use = 1:n.pcs,
    resolution = res.used, print.output = 0, save.SNN = TRUE,force.recalc=T) #DONE


# ```
#
# To visualize
# ```{r}
# If cells are too spread out, you can raise the perplexity.
#If you have few cells, try a lower perplexity (but never less than 10).
# https://distill.pub/2016/misread-tsne/
tiss <- RunTSNE(object = tiss, dims.use = 1:n.pcs, seed.use = 10, perplexity=70,
                check_duplicates = F)
# set check_duplicates = F when testing small number of genes since some cells might have identical profiles

# ```
#
# ```{r}
TSNEPlot(tiss, group.by = 'tissue')
# ```
#
# ```{r}
save(tiss, file=here("00_data_ingest", "11_global_robj", "FACS_all_preannotation.Robj"))
#load(here("00_data_ingest", "11_global_robj", "FACS_all_preannotation.Robj"))
```

# Add in metadata, annotations, and save.

# plot the location of a list of genes in the main tSNE plot
x11();FeaturePlot(tiss,bmp.ligands, cols.use = c("lightgrey","blue"),nCol=4  )

tiss@meta.data['cluster'] <- tiss@ident
tiss@meta.data['cell'] <- rownames(tiss@meta.data)

#anno = read_csv(here("00_data_ingest", "18_global_annotation_csv", "annotations_FACS.csv"))
anno = read_csv("/home/agranado/MEGA/Caltech/rnaseq/tabula-muris/00_data_ingest/18_global_annotation_csv/annotations_facs.csv")
anno %>% group_by(tissue, cell_ontology_class) %>% summarize(count = n())

# # A tibble: 120 x 3
# # Groups:   tissue [?]
#    tissue            cell_ontology_class                  count
#    <chr>             <chr>                                <int>
#  1 Aorta             endothelial cell                       188
#  2 Aorta             erythrocyte                             91
#  3 Aorta             fibroblast                              70
#  4 Aorta             professional antigen presenting cell    59
#  5 Bladder           bladder cell                           695
#  6 Bladder           bladder urothelial cell                683
#  7 Brain_Myeloid     macrophage                              61

tissue_colors = read_csv("/home/agranado/MEGA/Caltech/rnaseq/tabula-muris/00_data_ingest/15_color_palette/tissue_colors.csv")
colnames(tissue_colors) = c("tissue", "color")
#update metadata
# left_join()
# return all rows from x, and all columns from x and y. Rows in x with no match in y will have NA values in the new columns.
# If there are multiple matches between x and y, all combinations of the matches are returned.
tiss@meta.data <- tiss@meta.data %>%
		   left_join(anno %>% select(cell_ontology_class,cell_ontology_id,free_annotation, cell), by = 'cell') %>%
		   left_join(tissue_colors, by = 'tissue')

rownames(tiss@meta.data) <- tiss@meta.data$cell





#Here we can analyze the composition of each cluster:

#----let's take the important variables:
tiss@meta.data %>% select(cluster,cell_ontology_class,tissue) -> clusters.info
# > head(clusters.info)
#                       cluster cell_ontology_class        tissue
# A1.B000167.3_56_F.1.1       9          basal cell Mammary_Gland
# A1.B000412.3_56_F.1.1       3    endothelial cell         Heart
# A1.B000610.3_56_F.1.1      10        bladder cell       Bladder
# A1.B000633.3_56_F.1.1      13           leukocyte         Heart
# A1.B000826.3_39_F.1.1      11     microglial cell Brain_Myeloid
# A1.B001176.3_56_F.1.1       8     microglial cell Brain_Myeloid

#Lots of NA values in the cell_ontology_class : make it a string
 clusters.info$cell_ontology_class[is.na(clusters.info$cell_ontology_class)]<-"NA"


 #---- Order with respect to the count column
 cluster.size<-clusters.info %>% group_by(cell_ontology_class) %>% summarise(count=n())
 cluster.size = cluster.size[order(cluster.size$count,decreasing=F),]
 cluster.size$cell_ontology_class<-factor(cluster.size$cell_ontology_class,levels=cluster.size$cell_ontology_class)
 #barplot with the size of each cluster:
 p = ggplot(data = cluster.size,aes(x=cell_ontology_class,y=log10(count))) + geom_bar(stat="identity")
 x11();p + coord_flip()
#some cell types are represented by <100 cell, we can remove those

 cluster.size$freq <-cluster.size$count/sum(cluster.size$count)
 cluster.size = cluster.size[order(cluster.size$freq,decreasing=T),]


#-----we count how many cells are on each cluster and the relative frequency per cluster:
clusters.info %>% group_by(cluster,cell_ontology_class) %>% #Count how many cells have same cluster AND ontology
    summarise(count =n()) %>% group_by(cluster) %>% #summarise ungroups the data, so we group it again
        mutate(freq = count / sum(count)) -> clusters.count


      #   cluster cell_ontology_class     count    freq
      #   <fct>   <chr>                   <int>   <dbl>
      # 1 0       astrocyte                   8 0.00290
      # 2 0       B cell                     74 0.0268
      # 3 0       basal cell                219 0.0793
      # ....


#-----let's look at a particular cluster
      #What fraction of the cluster is a given cell type?
      clus.number = 10
      clus28<-clusters.count[clusters.count$cluster==clus.number,]
      clus28<-clus28[order(clus28$freq,decreasing=T),]
      clus28
    #   cluster cell_ontology_class              count   freq
    #   <fct>   <chr>                            <int>  <dbl>
    # 1 28      microglial cell                     62 0.223
    # 2 28      fibroblast                          49 0.176
    # 3 28      endothelial cell                    32 0.115
    # 4 28      mesenchymal cell                    19 0.0683
    # 5 28      mesenchymal stem cell of adipose    18 0.0647
    # 6 28      leukocyte                           12 0.0432
    # 7 28      stromal cell                        11 0.0396
    # 8 28      bladder cell                        10 0.0360

#-----A similar question: What fraction of a given cell type is in a particular cluster?
clusters.info %>% group_by(cluster,cell_ontology_class) %>% #Count how many cells have same cluster AND ontology
    summarise(count =n()) %>% group_by(cell_ontology_class) %>% #summarise ungroups the data, so we group it again
        mutate(freq = count / sum(count)) -> ontology.count

#--- Order by frequency
 cell.type = "microglial cell"
 astro.count = ontology.count[ontology.count$cell_ontology_class==cell.type,]
 astro.count=astro.count[order(astro.count$freq,decreasing=T),]
 astro.count$cluster<-factor(astro.count$cluster,levels=astro.count$cluster) #for gplot to respect the order
 head(astro.count)
 #
 # 1 21      astrocyte             141 0.675
 # 2 2       astrocyte              12 0.0574
 # 3 0       astrocyte               8 0.0383
 # 4 1       astrocyte               7 0.0335
 # 5 23      astrocyte               7 0.0335
 # 6 6       astrocyte               6 0.0287

#---- Look at how a given cell type is distributed across clusters:
 clusters.count[clusters.count$cell_ontology_class=="astrocyte",]




#---- Let's plot with the names of the overrepresented cell types
FetchData(tiss,c('tSNE_1','tSNE_2')) %>%
  summarize(xmin=min(tSNE_1),xmax=max(tSNE_1),ymin=min(tSNE_2),ymax=max(tSNE_2) )

  plot_min = min(lims$xmin, lims$ymin)
  plot_max = max(lims$xmax, lims$ymax)



tiss_FACS = tiss
save(tiss_FACS, file=here("00_data_ingest", "11_global_robj", "FACS_all.Robj"))
```
