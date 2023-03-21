source("./scripts/analysis/imports_new.R")

all_csv <- read.csv("./data/raw_data/pathbank/pathbank_all_proteins.csv")

pathway_filter <- all_csv %>% 
  filter(Species == "Mus musculus") %>%
  filter((Pathway.Subject == "Protein") | 
           (Pathway.Subject == "Signaling")) %>%
  group_by(Pathway.Name) %>%
  count() %>%
  filter(n>5) %>% 
  select(Pathway.Name) %>% 
  unique() %>%
  pull()

all_csv <- all_csv %>% 
  filter(Pathway.Name %in% pathway_filter) %>%
  select(Pathway.Name, Gene.Name)

colnames(all_csv) <- c("pathway", "gene")

our_pathways <- read.csv("./data/raw_data/pathbank/our_pathways.csv",
                         row.names = 1)

pathway_df <- rbind(all_csv, our_pathways)

pathway_df$pathway[pathway_df$pathway=='Notch'] <- 'Notch receptors, Dll ligands, and Fringe proteins'
pathway_df$pathway[pathway_df$pathway=='Bmp_Tgfb'] <- 'Tgf-beta family receptors'
pathway_df$pathway[pathway_df$pathway=='Srsf'] <- 'SRSF Splicing Protein Family'
pathway_df$pathway[pathway_df$pathway=='Eph_r'] <- 'Eph A-B receptors'
pathway_df$pathway[pathway_df$pathway=='Eph_l'] <- 'Ephrins'
pathway_df$pathway[pathway_df$pathway=='Wnt'] <- 'Frizzled and Lrp5 6 receptors for Wnt B Catenin Signaling'
pathway_df$pathway[pathway_df$pathway=='Fgfr'] <- 'FGF signaling'

write.csv(pathway_df, "./data/raw_data/pathbank/pathway_df.csv")
