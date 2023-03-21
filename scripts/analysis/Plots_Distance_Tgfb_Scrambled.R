set.seed(17) # Set seed to keep scrambling of the dataframe the same.
source("./scripts/analysis/imports_new.R")

## Change these variables
filename = "tgfb_trial" # how name is appended to files for pathway

## Directories for the files
output_dir = paste("./scripts/analysis/outputs/", filename, "_analysis/", filename, sep = "")
fig_dir = "./scripts/figures/"

pathway_df = read.csv('./data/raw_data/pathbank/pathway_df.csv', row.names = 1)

pathway_name =  'Tgf-beta family receptors' # tgfb pathway name
pathway_genes = genesPathway(pathway_name = pathway_name,
                             pathway_df = pathway_df,
                             seurat_obj = master_seurat)

min_genes_pathway = 2 # tgfb min. number of genes expressed
min_expr_threshold = 0.3 # tgfb minimum expression threshold for gene to be on
diverse_quantile = 0.9
optimal_k_pathway = 30 # Can be anything

## Returns MinMax normalized counts along with metadata at specified saturating quantile.
normalizedDevel <- function(seurat_obj = c(),    # Which Seurat object to use
                            pathway_genes = c(),           # List of pathway genes
                            sat_val = 0.99,                # Saturating quantile
                            fill_zero_rows = F             # If a gene has all 0's, fill with very small number
){
  
  # Get the log norm gene counts and annotations
  data_frame <- makeMainDataFrame(pathway_genes, 
                                  seurat_obj)
  
  if(!'cell_id' %in% names(seurat_obj@meta.data))
    data_frame %>% mutate(cell_id = paste(global_cluster, 
                                          dataset, 
                                          sep = "_")
    ) -> data_frame
  
  # Get dataframe of only gene expression values in the subset of cell states we want
  x = data_frame[, pathway_genes]
  
  # Compute the saturating values at a specified sat_val quantile
  max_sat_gene = apply(x, 2, quantile, sat_val) # starts from x
  
  # If there are genes that saturate to a value of 0, then replace the saturating value with the max. value of its expression
  if(sum(max_sat_gene == 0) > 0){
    max_val_gene = apply(x, 2, max) # starts from x
    max_sat_gene[max_sat_gene == 0] <- max_val_gene[max_sat_gene == 0]
  }
  
  # Actually saturate the values now
  for(s in 1:dim(x)[2])
    x[which(x[, s] > max_sat_gene[s]), s]<- max_sat_gene[s]
  
  # Apply MinMax Scaling to scale values between 0 to 1
  x <- minMaxNorm(x)
  
  # We can't do cosine distance on all zero genes so we just fill these with a very small value if we want
  if(fill_zero_rows)
    x[x == 0] = 10^-10
  
  # Modify the dataframe with the metadata to now have the MinMax scaled, saturated counts
  data_frame[, pathway_genes] <- x
  
  data_frame[,pathway_genes] <- randomizeColumns(df = data_frame[,pathway_genes],
                                                 pathway_genes=pathway_genes)
  
  row.names(data_frame) <- data_frame$global_cluster # Reset column names
  
  return(data_frame)
}

## Make a plot of k_opt profiles for cell types with pathway "ON
pipe_run <- quickPipeline(seurat_obj = master_seurat,
                          pathway_genes = pathway_genes,
                          k_final = optimal_k_pathway, 
                          min_genes_on = min_genes_pathway, 
                          min_expr = min_expr_threshold
)