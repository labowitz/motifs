source("./scripts/analysis/imports.R")
output_dir <- "./scripts/figures/figure_s1c/"

sce <- readRDS("./data/processed_data/forelimb.RDS")

# condition is gene1 
cond_two_genes = function(gene1, gene2, counts =c(), rand = F ){
  
  two_genes = counts[c(gene1,gene2),]
  if(rand){
    two_genes[1,] = sample(two_genes[1,])
    two_genes[2,] = sample(two_genes[2,])
  }
  pr_both = sum(two_genes[1,]>0 & two_genes[2,]>0)/dim(two_genes)[2]
  pr_gene1 = sum(two_genes[1,]>0)/dim(two_genes)[2]
  
  return(pr_cond = pr_both * pr_gene1)
  
}

tgfb_genes <- c("Acvr1", "Acvrl1", "Bmpr1a", "Bmpr1b", "Acvr2a", "Acvr2b", "Bmpr2")

# First cell type
cell_type = "Epithelial 2 "
filename = "epithelial_2"
pathway_name = "Tgfb"

assay(sce[tgfb_genes,which(colData(sce)$cell_type == cell_type)], "X") -> counts_2

p_mat = matrix(0, length(tgfb_genes), length(tgfb_genes))
x = counts_2 

# condition is in the rows P(g2 | g1) 
for(i in 1:length(tgfb_genes)){
  gene1 = tgfb_genes[i]
  for(j in 1:length(tgfb_genes)){
    if(i != j){ # only if the genes are different 
      gene2 = tgfb_genes[j]
      # repeats bootstrap 
      rand_cond = rep(0,1000)
      
      for(r in 1:1000)
        rand_cond[r] = cond_two_genes(gene1, gene2, x,rand=T)
      # real pr for this pair
      p_real = cond_two_genes(gene1, gene2, x,rand=F)
      # p_value based on the distribution 
      p_val = sum(rand_cond >= p_real)/length(rand_cond)
      
      if(p_val==0) p_val = 0.0001 # as limited by the sample size
      p_mat[i,j] = -log10(p_val)
    }
  }
}

row.names(p_mat) <- tgfb_genes
colnames(p_mat) <- tgfb_genes

# annotation for significant pairs 
text_ann = matrix('',dim(p_mat)[1],dim(p_mat)[2])
text_ann[p_mat>=-log10(0.05)] = '*'
text_ann[p_mat>=-log10(0.001)] = '**'

# heatmap 
x_mat = p_mat
x_mat[x_mat==0] = NA

pdf(paste(output_dir, filename, "_", pathway_name,".pdf", sep = ""))
superheat(x_mat, 
          scale = F, 
          X.text = text_ann ,
          heat.na.col = 'gray',
          pretty.order.rows = F, 
          pretty.order.cols = F, 
          title='Epithelial 2')
dev.off()

ct_rows = master_seurat@meta.data %>% filter(grepl("Forelimb", dataset) & grepl("Epithelial 2", cell_ontology_class)) %>% rownames()
df <- normalizedDevel(this_pathway = genesPathway("Bmp_Tgfb"), master_seurat = master_seurat)
df <- df %>% filter(cell_id %in% ct_rows) %>% select(tgfb_genes)

df <- data.frame(counts = colMeans(df), names = colMeans(df) %>% names())
p<-ggplot(data=df, aes(x = names, y = counts)) +
  geom_bar(stat="identity")

p <- p + labs(x = "gene", y = "Norm. Counts") + ggtitle("Epithelial 2") + 
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 25)),
               axis.text.x = element_text(angle = 90, size = 18, vjust = 0.5, margin = margin(t = -10)), 
               axis.ticks.x = element_blank(),
               axis.title.y = element_text(size = 18, margin = margin(r = 25)),
               axis.text.y = element_text(size = 18),
               axis.ticks.y = element_blank(),
               panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               panel.background = element_blank(),
               plot.title = element_text(size = 20)
        )

ggsave(paste(output_dir, filename, "_", pathway_name,"_barplot.pdf", sep = ""))

# Second cell type
# First cell type
cell_type = "Mesenchymal 2 "
filename = "mesenchymal_2"
pathway_name = "Tgfb"

assay(sce[tgfb_genes,which((colData(sce)$cell_type == cell_type) & (colData(sce)$stage == 15.0))], "X") -> counts_2

p_mat = matrix(0, length(tgfb_genes), length(tgfb_genes))
x = counts_2 

# condition is in the rows P(g2 | g1) 
for(i in 1:length(tgfb_genes)){
  gene1 = tgfb_genes[i]
  for(j in 1:length(tgfb_genes)){
    if(i != j){ # only if the genes are different 
      gene2 = tgfb_genes[j]
      # repeats bootstrap 
      rand_cond = rep(0,1000)
      
      for(r in 1:1000)
        rand_cond[r] = cond_two_genes(gene1, gene2, x,rand=T)
      # real pr for this pair
      p_real = cond_two_genes(gene1, gene2, x,rand=F)
      # p_value based on the distribution 
      p_val = sum(rand_cond >= p_real)/length(rand_cond)
      
      if(p_val==0) p_val = 0.0001 # as limited by the sample size
      p_mat[i,j] = -log10(p_val)
    }
  }
}

row.names(p_mat) <- tgfb_genes
colnames(p_mat) <- tgfb_genes

# annotation for significant pairs 
text_ann = matrix('',dim(p_mat)[1],dim(p_mat)[2])
text_ann[p_mat>=-log10(0.05)] = '*'
text_ann[p_mat>=-log10(0.001)] = '**'

# heatmap 
x_mat = p_mat
x_mat[x_mat==0] = NA

pdf(paste(output_dir, filename, "_", pathway_name,".pdf", sep = ""))
superheat(x_mat, 
          scale = F, 
          X.text = text_ann ,
          heat.na.col = 'gray',
          pretty.order.rows = F, 
          pretty.order.cols = F, 
          title='Epithelial 2')
dev.off()

ct_rows = master_seurat@meta.data %>% filter(grepl("Forelimb", dataset) & grepl("Mesenchymal 2", cell_ontology_class) & grepl("E15", age)) %>% rownames()
df <- normalizedDevel(this_pathway = genesPathway("Bmp_Tgfb"), master_seurat = master_seurat)
df <- df %>% filter(cell_id %in% ct_rows) %>% select(tgfb_genes)

df <- data.frame(counts = colMeans(df), names = colMeans(df) %>% names())
p<-ggplot(data=df, aes(x = names, y = counts)) +
  geom_bar(stat="identity")

p <- p + labs(x = "gene", y = "Norm. Counts") + ggtitle("Mesenchymal 2") + 
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 25)),
        axis.text.x = element_text(angle = 90, size = 18, vjust = 0.5, margin = margin(t = -10)), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18, margin = margin(r = 25)),
        axis.text.y = element_text(size = 18),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 20)
  )

ggsave(paste(output_dir, filename, "_", pathway_name,"_barplot.pdf", sep = ""))

# Third cell type
cell_type = "Mesenchymal 2 "
filename = "mesenchymal_2"
pathway_name = "Tgfb"

assay(sce[tgfb_genes,which((colData(sce)$cell_type == cell_type) & (colData(sce)$stage == 15.0))], "X") -> counts_2

p_mat = matrix(0, length(tgfb_genes), length(tgfb_genes))
x = counts_2 

# condition is in the rows P(g2 | g1) 
for(i in 1:length(tgfb_genes)){
  gene1 = tgfb_genes[i]
  for(j in 1:length(tgfb_genes)){
    if(i != j){ # only if the genes are different 
      gene2 = tgfb_genes[j]
      # repeats bootstrap 
      rand_cond = rep(0,1000)
      
      for(r in 1:1000)
        rand_cond[r] = cond_two_genes(gene1, gene2, x,rand=T)
      # real pr for this pair
      p_real = cond_two_genes(gene1, gene2, x,rand=F)
      # p_value based on the distribution 
      p_val = sum(rand_cond >= p_real)/length(rand_cond)
      
      if(p_val==0) p_val = 0.0001 # as limited by the sample size
      p_mat[i,j] = -log10(p_val)
    }
  }
}

row.names(p_mat) <- tgfb_genes
colnames(p_mat) <- tgfb_genes

# annotation for significant pairs 
text_ann = matrix('',dim(p_mat)[1],dim(p_mat)[2])
text_ann[p_mat>=-log10(0.05)] = '*'
text_ann[p_mat>=-log10(0.001)] = '**'

# heatmap 
x_mat = p_mat
x_mat[x_mat==0] = NA

pdf(paste(output_dir, filename, "_", pathway_name,".pdf", sep = ""))
superheat(x_mat, 
          scale = F, 
          X.text = text_ann ,
          heat.na.col = 'gray',
          pretty.order.rows = F, 
          pretty.order.cols = F, 
          title='Epithelial 2')
dev.off()

ct_rows = master_seurat@meta.data %>% filter(grepl("Forelimb", dataset) & grepl("Mesenchymal 2", cell_ontology_class) & grepl("E15", age)) %>% rownames()
df <- normalizedDevel(this_pathway = genesPathway("Bmp_Tgfb"), master_seurat = master_seurat)
df <- df %>% filter(cell_id %in% ct_rows) %>% select(tgfb_genes)

df <- data.frame(counts = colMeans(df), names = colMeans(df) %>% names())
p<-ggplot(data=df, aes(x = names, y = counts)) +
  geom_bar(stat="identity")

p <- p + labs(x = "gene", y = "Norm. Counts") + ggtitle("Mesenchymal 2") + 
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 25)),
        axis.text.x = element_text(angle = 90, size = 18, vjust = 0.5, margin = margin(t = -10)), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18, margin = margin(r = 25)),
        axis.text.y = element_text(size = 18),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 20)
  )

ggsave(paste(output_dir, filename, "_", pathway_name,"_barplot.pdf", sep = ""))