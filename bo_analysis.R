# Imports
source("./scripts/analysis/imports.R")

## Change these variables
filename = "mediator" # how name is appended to files for pathway
gene_list <- c("Med1",
               "Med9",
               "Med19",
               "Med15",
               "Med16",
               "Med23",
               "Med24",
               "Med25",
               "Med27",
               "Med28",
               "Med29",
               "Med30")

pathway_name =  "Mediator" # Tgfb pathway name
min_genes_pathway = 2 # Tgfb min. number of genes expressed
min_expr_threshold = 0.2 # Tgfb minimum expression threshold for gene to be on

## Directories for the files
output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")
fig_dir = "./scripts/figures/"

## Add the gene list to our list of pathway genes we reference.
all_pathways <- rbind(all_pathways, data.frame(pathway = rep(pathway_name, length(gene_list)), gene = gene_list))

optimal_k_pathway = 30 # optimal pathway components, computed from z-score

## Make a plot of k_opt profiles for cell types with pathway "ON
pipe_run <- quickPipeline(master_seurat = master_seurat, 
                          k_final = optimal_k_pathway, 
                          silent_plot = F, 
                          which_pathway = pathway_name, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold, 
                          which_profiles = 'both', 
                          save_pdf = T, 
                          pdf_file = paste(fig_dir, "Mediator.pdf", sep = ""))

