library(shiny)
library(DIZtools)

setwd("~/Documents/Elowitz Lab/pathway_motifs")
source("module/helper.r")
options(dplyr.summarise.inform = FALSE)

## Seurat object to import
master_seurat <- readRDS("./data/processed_data/master_seurat.RDS") # Most recent Seurat object
pca_proj <- readRDS("./data/processed_data/umap_may.RDS") # PCA projection from 1928 HVGs
hvg_genes <- readRDS("./data/processed_data/hvg_genes.RDS") # list of HVGs (shouldn't be used)


## Set of files to import 
pathway_list <- readRDS("./data/processed_data/pathway_list.RDS") # List of pathway components
param_list <- readRDS("./data/processed_data/param_list.RDS") # List of pathway parameters
all_pathways <- readRDS("./data/processed_data/all_pathways.RDS") # List of all pathways


## Colors to import
colors_1206 <- readRDS("./data/processed_data/colors_1206.RDS") # List with labels and corresponding colors
glasbey <- readRDS("./data/processed_data/glasbey.RDS") # List of glasbey colors for categorical labels


# Alter the metadata of the file to reflect changes in "Organ_specific" labels
organ_specific <- read.table('./data/raw_data/organ_specific_metadata - organ_specific_metadata.tsv', header = T, sep = "\t")

# Save previous 
master_seurat@meta.data$Cell_class_prev <- master_seurat@meta.data$Cell_class

# Check indexes are the same 
#identical(row.names(master_seurat@meta.data[ master_seurat@meta.data$Cell_class=='Organ_specific',  ]) , organ_specific$index ) 

# replace with new Cell type classes 
master_seurat@meta.data[master_seurat@meta.data$Cell_class == 'Organ_specific', ]$Cell_class <- organ_specific$Cell_class

ui <- fluidPage(
  
  # Application title
  titlePanel("Pathway motif Browser"),
  
  fluidRow(
    
    column(4, 
           helpText("Enter pathway name."),
           textInput(
             inputId = "pathway_name",
             label = "Enter pathway name", 
             value = "Bmp_Tgfb"),
           helpText("Enter the pathway name and gene list."),
           selectizeInput(
             inputId = "geneList",
             label = "Enter gene names",
             choices = row.names(master_seurat),
             selected = all_pathways[all_pathways$pathway=="Bmp_Tgfb",]$gene,
             multiple = TRUE,
             width = "100%",
             options = list(
               'plugins' = list('remove_button'),
               'create' = TRUE,
               'persist' = TRUE
             )
           ),
           helpText("Enter the minimum number of genes for the pathway to be considered 'ON'."),
           numericInput(
             inputId = "min_genes_pathway",
             label = "Enter min. # of genes on",
             value = 2,
           ),
           helpText("Enter the minimum expression value for a gene to be considered 'ON' from a scale of 0 to 0.1"),
           numericInput(
             inputId = "min_expr_threshold",
             label = "Enter min. expression threshold",
             min = 0.0,
             max = 1.0,
             value = 0.2,
           ),
           helpText("Press when you are ready to perform the data analysis!"),
           actionButton(
             inputId = "ready",
             label = "Press when ready!",
           ),
    ),
    column(4,
           mainPanel(
             plotOutput("silh_plot"),
             downloadButton("download_silh_plot",
                            label="Download Silhouette Z-Score Plot")
           )
    ),
    column(4, 
           mainPanel(
             plotOutput("heatmap")
           )
    )
  ),
  fluidRow(
    column(6,
           mainPanel(
             plotOutput("ecdf_plot")
           )
    ),
    column(6,
           mainPanel(
             plotOutput("rank_plot")
           )
    ),
  ),
  fluidRow(
    column(6,
           mainPanel(
             plotOutput("motif_heatmap")
           )
    ),
    column(6,
           mainPanel(
             plotOutput("motif_type")
           )
    ),
  )
)

server <- function(input, output, session) {
  
  v <- reactiveValues(
    silh_plt = readRDS("./data/processed_data/tgfb_silh_plt.RDS"),
    silh_z = readRDS("./data/processed_data/tgfb_silh_z.RDS")[[1]],
    control_res = readRDS("./data/processed_data/tgfb_control_res.RDS"),
    opt_k = readRDS("./data/processed_data/tgfb_silh_z.RDS")[[2]],
  )
  
  observeEvent(
    input$ready,
    {
      v$silh_plt = silhouettePlot(which_pathway = input$pathway_name, 
                                  min_ON = input$min_genes_pathway, 
                                  min_expr = input$min_expr_threshold, 
                                  n_bootstraps = 10
      )
      temp = silhouette_zscore(silh_result = v$silh_plt,
                                   pathway_name = input$pathway_name,
                                   x_offset = 6,
                                   max.y = 0.6,
                                   min.y = 0.1 # adjust axis parameters 
      )
      
      v$silh_z = temp[[1]]
      
      v$opt_k = temp[[2]]
      
      v$control_res = fullControlPathway(this_pathway = input$pathway_name,
                                         k_pathway = v$opt_k , 
                                         filter_pathway = 'both', # adult + devel 
                                         this_seurat = master_seurat, # seurat object
                                         null_list = hvg_genes, #list of highly variable genes 
                                         n_samples = 100, 
                                         filter_manual = T,
                                         min_expr_gene = input$min_expr_threshold, 
                                         min_genes_ON = input$min_genes_pathway, # default > 2 genes ON 
                                         n_pcs_global = 100, # how many PCs to use
                                         embedding_matrix = pca_proj # PCA embedding for diversity 
      )
    }
  )
  
  output$silh_plot <- renderPlot({
    
    v$silh_z
    
  })
  
  output$heatmap <- renderPlot({
    
    quickPipeline(
      master_seurat = master_seurat, 
      k_final = v$opt_k, 
      which_pathway = input$pathway_name, 
      min_genes_on = input$min_genes_pathway, 
      min_expr = input$min_expr_threshold, 
      which_profiles = 'both'
    )[[3]]
    
  })
  
  output$ecdf_plot <- renderPlot({
    ecdf_diversity(v$control_res,
                   pathway_name=input$pathway_name)
  })
  
  output$rank_plot <- renderPlot({
    rank_diversity(
      which_pathway = input$pathway_name,
      k_motifs = v$opt_k,
      min_expression = input$min_expr_threshold,
      min_genes_pathway = input$min_genes_pathway,
      embedding_matrix = pca_proj
    )
  })
  
  output$motif_heatmap <- renderPlot({
    ## Plot the most diverse profiles, aka the motifs
    # 1. New: Select all profiles 
    diverse_df <- v$control_res$profiles
    # 2. tidy data.frame to average gene expression within a pathway profile 
    diverse_df %>% pivot_longer(cols = genesPathway(input$pathway_name), 
                                names_to = 'gene', 
                                values_to = 'expression') %>% 
      select(cell_id, cell_ontology_class, Tissue, Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy 
    
    # 3. average profile 
    diverse_tidy %>% group_by(class_label,gene) %>% summarise(mean_expr = mean(expression), 
                                                              rank = mean(rank),diversity = mean(diversity), 
                                                              cell_types = mean(n)) -> diverse_tidy
    
    # 4. wide format 
    diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types), names_from=gene,values_from=mean_expr) %>% tibble::column_to_rownames('class_label')-> diverse_mat
    
    # 5. We do the filtering here either for motifs or for non-diverse profiles 
    v$control_res$diversity %>% 
      dplyr::filter(type=='transcriptome') %>% # choose the null model 
      pull(d) %>% quantile(0.90) -> divers_quantile
    
    diverse_mat %>% dplyr::filter(diversity>divers_quantile) -> diverse_mat 
    
    #pdf(paste(output_dir, "_motif_profiles.pdf", sep = ""))
    motif_heatmap <- superheat(diverse_mat[,genesPathway(input$pathway_name)],
                               pretty.order.rows = T,
                               heat.pal = black_pal(10),
                               bottom.label.text.angle = 90, 
                               yr = sqrt(diverse_mat$cell_types),
                               yr.plot.type='bar',
                               yr.axis.name = "N cell types",
                               row.title = "Pathway motifs",
                               column.title = "Pathway genes",
                               bottom.label.size = 0.3,
                               grid.hline.col = "white",
                               grid.hline.size = 2, 
                               
    )
    #dev.off()
  })
  
  output$motif_type <- renderPlot({
    ## Plot the most diverse profiles, aka the motifs
    # 1. New: Select all profiles 
    diverse_df <- v$control_res$profiles
    # 2. tidy data.frame to average gene expression within a pathway profile 
    diverse_df %>% pivot_longer(cols = genesPathway(input$pathway_name), 
                                names_to = 'gene', 
                                values_to = 'expression') %>% 
      select(cell_id, cell_ontology_class, Tissue, Cell_class, gene, expression, class_label, rank, diversity, n) -> diverse_tidy 
    
    # 3. average profile 
    diverse_tidy %>% group_by(class_label,gene) %>% summarise(mean_expr = mean(expression), 
                                                              rank = mean(rank),diversity = mean(diversity), 
                                                              cell_types = mean(n)) -> diverse_tidy
    
    # 4. wide format 
    diverse_tidy %>% pivot_wider(id_cols = c(class_label,rank, diversity, cell_types), names_from=gene,values_from=mean_expr) %>% tibble::column_to_rownames('class_label')-> diverse_mat
    
    # 5. We do the filtering here either for motifs or for non-diverse profiles 
    v$control_res$diversity %>% 
      dplyr::filter(type=='transcriptome') %>% # choose the null model 
      pull(d) %>% quantile(0.90) -> divers_quantile
    
    diverse_mat %>% dplyr::filter(diversity>divers_quantile) -> diverse_mat 
    
    
    # Motif profile cell type distribution
    diverse_df %>% dplyr::filter(class_label %in% row.names(diverse_mat)) %>% 
      select(Tissue, class_label ) %>% dplyr::filter(!Tissue %in% c('mat','endoderm')) %>% group_by(class_label,Tissue)%>% 
      count  %>% 
      pivot_wider(id_cols=class_label,names_from = Tissue,
                  values_from= n,values_fill = 0) -> tissue_distribution
    
    
    tissue_pal<-colorRampPalette(brewer.pal(n = 9, name = 'PuRd'))
    x = tissue_distribution %>% ungroup %>% select(-class_label) %>% as.matrix() 
    row.names(x) <- tissue_distribution$class_label
    # make tissue names pretty 
    colnames(x)<-str_replace(string = colnames(x),pattern ="_",replacement = ' ') %>% DIZtools::firstup() 
    
    x <- x[,sort(colnames(x))]
    
    pheatmap(sqrt(x), treeheight_row = 20,treeheight_col = 20,
             clustering_method = 'ward.D2',col = tissue_pal(100),
             cluster_rows = F,
             cluster_cols = F,
             fontsize =12,angle_col = 45,
             height = 4, width = 6) # keep this size to make it similar to the motif heatmap 
    
  })
  
  output$download_silh_plot <- downloadHandler(
    filename = function() {
      paste0("silh_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file,
             output$silh_plot,
             width = 7,
             height = 5)
    }
  )
}

shinyApp(ui = ui, server = server)