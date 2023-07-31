library(shiny)

setwd("~/Documents/Elowitz Lab/pathway_motifs")
source("module/helper.R")
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
             plotOutput("silh_plot")
           )
    ),
    column(4, 
           mainPanel(
             plotOutput("heatmap")
           )
    )
  ),
)

server <- function(input, output, session) {
  
  v <- reactiveValues(
    silh_plt = readRDS("data/processed_data/tgfb_silh_plt.RDS"),
    silh_z = readRDS("data/processed_data/tgfb_silh_z.RDS")[[1]],
    #control_res = readRDS("data/processed_data/tgfb_control_res.RDS"),
    opt_k = readRDS("data/processed_data/tgfb_silh_z.RDS")[[2]],
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
  
}

shinyApp(ui = ui, server = server)