library(dplyr)
library(ggplot2)
library(gridExtra)

# Load Pathbank analysis results
load("~/Documents/Elowitz Lab/pathway_motifs/data/processed_data/pathbank_silhouetteControl2022.rda")

# Pathway genes list
all_pathways <- read.table("./data/raw_data/allPathways_listGenes_dec2021.tsv", header = T, sep = "\t")

all_pathways$pathway[all_pathways$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
all_pathways$pathway[all_pathways$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
all_pathways$pathway[all_pathways$pathway=='Srsf'] <- 'SRSF Splice Regulators'
all_pathways$pathway[all_pathways$pathway=='Eph_r'] <- 'Eph-Ephrin'
all_pathways$pathway[all_pathways$pathway=='Eph_l'] <- 'Eph-Ephrin'
all_pathways$pathway[all_pathways$pathway=='Wnt'] <- 'Frizzled and Lrp5/6 receptors for Wnt B/Catenin Signaling'
all_pathways$pathway[all_pathways$pathway=='Fgfr'] <- 'FGF cell signaling proteins'

data_real_plot$pathway[data_real_plot$pathway=='Srsf'] <- 'SRSF Splice Regulators'
data_real_plot$pathway[data_real_plot$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
data_real_plot$pathway[data_real_plot$pathway=='Eph_r'] <- 'Eph-Ephrin'
data_real_plot$pathway[data_real_plot$pathway=='Eph_l'] <- 'Eph-Ephrin'
data_real_plot$pathway[data_real_plot$pathway=='Notch'] <- 'Notch receptors, Dll ligands and Fringe proteins'
data_real_plot$pathway[data_real_plot$pathway=='Wnt'] <- 'Frizzled and Lrp5/6 receptors for Wnt B/Catenin Signaling'
data_real_plot$pathway[data_real_plot$pathway=='Fgfr'] <- 'FGF cell signaling proteins'

# Plot densities
plotDensity <- function(){
  # Make master DF and plot 
  all_data_pathbank <- rbind(data_real_rank %>% dplyr::select(pathway, n_genes, pathway_class, diff, max_z), 
                             data_random_diff %>% dplyr::select(pathway, n_genes, pathway_class, diff, max_z)) 
  
  all_data_pathbank[all_data_pathbank$pathway_class!='random',]$pathway_class <- 'pathway'
  
  
  p5 <- all_data_pathbank %>% ggplot(aes(x = diff, fill = pathway_class)) +geom_density(alpha = 0.6) + theme_minimal()  + theme(legend.position = "none") + 
    geom_vline(xintercept = quantile(data_random_diff$diff, 0.01)) + 
    coord_cartesian(xlim = c(-4,5))
  
  p6 <- data_real_rank %>% ggplot(aes(x = diff)) + geom_histogram() + theme_minimal() +
    geom_vline(xintercept = quantile(data_random_diff$diff, 0.01)) + 
    geom_vline(xintercept = quantile(data_random_diff$diff, 0.99)) + 
    coord_cartesian(xlim = c(-4,5))
  return(list(p5,p6))	
}

# Individual histograms
clusterScoreDist <- function(which_pathway = 'Bmp_Tgfb'){
  pathway_size <- all_pathways %>% dplyr::filter(pathway== which_pathway) %>% pull(gene) %>% length 
  null_dist <-  data_random_plot %>% dplyr::filter(n_genes == pathway_size)
  
  p1 <- null_dist %>% ggplot(aes(x = cluster_ratio)) + 
    #geom_histogram(aes(y = ..density..)) + 
    geom_density(col = "blue") + 
    theme_minimal()
  
  path_stats <- data_real_plot %>% dplyr::filter(pathway==which_pathway) 
  
  p1 <- p1 + geom_vline(xintercept = path_stats$max_peak/path_stats$n_genes, col = 'blue') + xlab('N clusters / n genes')+ 
    theme(text = element_text(size = 10)) + 
    ggtitle(paste(which_pathway , '\nvs ', dim(null_dist)[1],' random sets' , sep =""))
  return(p1) 
}

pathway_order <- c(
  "CXCR4 Signaling Pathway",
  "Lysophosphatidic Acid LPA6 Signalling",
  "Rac 1 Cell Motility Signaling Pathway",
  "Cadmium Induces DNA Synthesis and Proliferation in Macrophages ",
  "Apoptotic DNA Fragmentation and Tissue Homeostasis",
  "GnRH Signaling Pathway",
  "Growth Hormone Signaling Pathway ",
  "SRSF Splice Regulators",
  "Notch receptors, Dll ligands and Fringe proteins",
  'Frizzled and Lrp5/6 receptors for Wnt B/Catenin Signaling',
  "Ubiquitin–Proteasome Pathway",
  'Tgf-beta family receptors',
  'Eph-Ephrin'
)

p1 <- clusterScoreDist(pathway_order[[1]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) + theme(axis.line = element_line(size = 0.5),
                                                                                                                        panel.grid.minor.x = element_blank(),
                                                                                                                        panel.grid.minor.y = element_blank())
p2 <- clusterScoreDist(pathway_order[[2]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p3 <- clusterScoreDist(pathway_order[[3]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p4 <- clusterScoreDist(pathway_order[[4]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p5 <- clusterScoreDist(pathway_order[[5]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p6 <- clusterScoreDist(pathway_order[[6]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p7 <- clusterScoreDist(pathway_order[[7]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p8 <- clusterScoreDist(pathway_order[[8]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p9 <- clusterScoreDist(pathway_order[[9]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                       axis.ticks.y = element_blank(),
                                                                                                                       axis.title.y = element_blank(),
                                                                                                                       axis.line = element_line(size = 0.3),
                                                                                                                       panel.grid.minor.x = element_blank(),
                                                                                                                       panel.grid.minor.y = element_blank())
p10 <- clusterScoreDist(pathway_order[[10]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                         axis.ticks.y = element_blank(),
                                                                                                                         axis.title.y = element_blank(),
                                                                                                                         axis.line = element_line(size = 0.3),
                                                                                                                         panel.grid.minor.x = element_blank(),
                                                                                                                         panel.grid.minor.y = element_blank())
p11 <- clusterScoreDist(pathway_order[[11]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                         axis.ticks.y = element_blank(),
                                                                                                                         axis.title.y = element_blank(),
                                                                                                                         axis.line = element_line(size = 0.3),
                                                                                                                         panel.grid.minor.x = element_blank(),
                                                                                                                         panel.grid.minor.y = element_blank())
p12 <- clusterScoreDist(pathway_order[[12]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                         axis.ticks.y = element_blank(),
                                                                                                                         axis.title.y = element_blank(),
                                                                                                                         axis.line = element_line(size = 0.3),
                                                                                                                         panel.grid.minor.x = element_blank(),
                                                                                                                         panel.grid.minor.y = element_blank())
p13 <- clusterScoreDist(pathway_order[[13]]) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+ theme(axis.text.y = element_blank(),
                                                                                                                         axis.ticks.y = element_blank(),
                                                                                                                         axis.title.y = element_blank(),
                                                                                                                         axis.line = element_line(size = 0.3),
                                                                                                                         panel.grid.minor.x = element_blank(),
                                                                                                                         panel.grid.minor.y = element_blank())

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, nrow = 1) -> g

ggsave("./scripts/figures/Figure_5_figure_supplement_1B.pdf", g, width = 50, height = 5, limitsize = F)