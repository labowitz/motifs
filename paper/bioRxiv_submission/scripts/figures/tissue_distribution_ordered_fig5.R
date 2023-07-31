tgfb_tissue <- read.csv("./scripts/figures/tgfb_tissue.csv")
wnt_tissue <- read.csv("./scripts/figures/wnt_tissue.csv")
notch_tissue <- read.csv("./scripts/figures/notch_tissue.csv")
eph_tissue <- read.csv("./scripts/figures/eph-ephrin_tissue.csv")

df_arr <- list(tgfb_tissue, wnt_tissue, notch_tissue, eph_tissue)

all_tissues <- union(union(union(colnames(tgfb_tissue), colnames(wnt_tissue)), colnames(notch_tissue)), colnames(eph_tissue))

for (i in sort(all_tissues[2:length(all_tissues)])) {
  for (j in 1:length(df_arr)) {
    if (!(i %in% colnames(df_arr[[j]]))) {
      df_arr[[j]][i] <- rep(0, dim(df_arr[[j]])[1])
    }
  }
}

tgfb_order <- c(21, 10, 24, 9, 8, 14, 13, 16, 7, 20, 15, 23, 22, 27)
wnt_order <- c(25, 10, 15, 16, 14, 4, 3, 28, 17, 12, 20, 26)
notch_order <- c(29, 7, 20, 9, 5, 23, 19, 16, 6, 4, 2, 30, 28, 12, 26, 17)
eph_order <- c(24, 18, 53, 20, 35, 26, 47, 23, 46, 38, 24, 52, 25, 51, 54, 37, 50, 49)

order_arr <- list(tgfb_order, wnt_order, notch_order, eph_order)

pathway_names <- c("tgfb", "wnt", "notch", "eph")

for (i in 1:length(df_arr)){
  x <- df_arr[[i]][match(order_arr[[i]], df_arr[[i]]$X),]
  x <- subset(x, select = -X)
  tissue_pal<-colorRampPalette(brewer.pal(n = 9, name = 'PuRd'))
  
  pheatmap(sqrt(x), treeheight_row = 20,treeheight_col = 20,
           clustering_method = 'ward.D2',col = tissue_pal(100),
           cluster_rows = F,
           cluster_cols = F,
           fontsize =12,angle_col = 45,
           filename = paste("./scripts/figures/", pathway_names[[i]], "_tissue_distribution.pdf", sep = ""),
           height =4, width = 6) # keep this size to make it similar to the motif heatmap 
  
}

