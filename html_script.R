library(shiny)
library(Cairo)

source("./html_module_new.R")
options(dplyr.summarise.inform = FALSE)

## Seurat object to import
master_seurat <- readRDS("./data/processed_data/master_seurat.RDS") # Most recent Seurat object
master_seurat@meta.data$Tissue[master_seurat@meta.data$Tissue=="intestine"] = "dev. intestine"
annotations <- createAnnotations("./data/processed_data/integrated_meta_data.csv")
all_pathways = read.csv("./data/raw_data/pathbank/pathway_df.csv", row.names = 1)

## Colors to import
colors_1206 <- readRDS("./data/processed_data/colors_1206.RDS") # List with labels and corresponding colors
glasbey <- readRDS("./data/processed_data/glasbey.RDS") # List of glasbey colors for categorical labels

ui <- fluidPage(
  
  # Application title
  titlePanel("Pathway Motif Browser"),
  
  fluidRow(
    
    column(3, 
           helpText("Enter pathway name."),
           textInput(
             inputId = "pathway_name",
             label = "Enter pathway name, no spaces please—use underscores!", 
             value = 'Tgf-beta family receptors'),
           helpText("Enter the pathway name and gene list."),
           selectizeInput(
             inputId = "geneList",
             label = "Enter gene names",
             choices = row.names(master_seurat),
             selected = all_pathways[all_pathways$pathway=='Tgf-beta family receptors',]$gene,
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
    column(3,
           mainPanel(
             plotOutput("silh_z"),
             downloadButton('d_silh_z', 'Download silhouette plot.')
           ),
           mainPanel(
             plotOutput("heatmap"),
             downloadButton('d_heatmap', 'Download heatmap plot.'),
             downloadButton('d_heatmap_data', 'Download heatmap data')
           )
    ),
    column(6, 
           mainPanel(
             plotOutput("glob_dendr_diversity"),
             downloadButton('d_glob_dendr_diversity', 'Download global heatmap plot.')
           ),
           mainPanel(
             plotOutput("global_umap"),
             downloadButton('d_global_umap', 'Download heatmap plot.')
           )
    )
  ),
  
  fluidRow(
    column(6,
           mainPanel(
             plotOutput("ecdf_diversity"),
             downloadButton('d_ecdf_diversity', 'Download diversity plot.')
           )
    ),
    column(6,
           mainPanel(
             plotOutput("rank_diversity"),
             downloadButton('d_rank_diversity', 'Download rank diversity plot.')
           )
    )
  ),
  
  fluidRow(
    column(6,
           mainPanel(
             plotOutput("motif_heatmap"),
             downloadButton('d_motif_heatmap', 'Download motif heatmap.')
           )
    ),
    
    column(6,
           mainPanel(
             plotOutput("motif_ct_heatmap"),
             downloadButton('d_motif_ct_heatmap', 'Download motif cell type heatmap.')
           )
    )
    
  )
)

server <- function(input, output, session) {
  
  v <- reactiveValues(all_pathways = read.csv("./data/raw_data/pathbank/pathway_df.csv", row.names = 1), # List of all pathways
                      pathway_genes = genesPathway(pathway_name = 'Tgf-beta family receptors',
                                                   pathway_df = read.csv("./data/raw_data/pathbank/pathway_df.csv"),
                                                   seurat_obj=master_seurat),
                      silh_plt = readRDS("./data/processed_data/tgfb_silh_plt.RDS"),
                      silh_z = readRDS("./data/processed_data/tgfb_silh_z.RDS")[[1]],
                      opt_k = readRDS("./data/processed_data/tgfb_silh_z.RDS")[[2]],
                      control_res = readRDS("./data/processed_data/tgfb_control_res.RDS")
  )
  
  observeEvent(
    input$ready,
    {
      v$all_pathways <- addGenes(pathway_df = v$all_pathways,
                                 pathway_genes = input$geneList, 
                                 pathway_name = input$pathway_name
      )
      
      v$pathway_genes <- genesPathway(pathway_name = input$pathway_name,
                                      pathway_df = v$all_pathways
      )
      
      v$silh_plt = silhouettePlot(pathway_genes = v$pathway_genes, 
                                  pathway_name = input$pathway_name,
                                  min_genes_on = input$min_genes_pathway, 
                                  min_expr = input$min_expr_threshold, 
                                  n_bootstraps = 10,
                                  seurat_obj = master_seurat
      )
      
      v$silh_z = silhouette_zscore(v$silh_plt,
                                   min_expr = input$min_expr_threshold,
                                   pathway_name = input$pathway_name)[[1]]
      
      v$opt_k = perc_k_finder(silhouette_zscore(v$silh_plt,
                              min_expr = input$min_expr_threshold,
                              pathway_name = input$pathway_name)[[2]],
                              percentile = 0.9)
      
      v$control_res = fullControlPathway(pathway_genes = v$pathway_genes,
                                         k_final = v$opt_k,
                                         seurat_obj = master_seurat, # seurat object
                                         null_list = hvg_genes, #list of highly variable genes 
                                         n_samples = 10, 
                                         filter_manual = T,
                                         min_genes_on = input$min_genes_pathway, 
                                         min_expr = input$min_expr_threshold, 
                                         n_pcs = 100, # how many PCs to use
                                         manual_embedding = pca_proj, # PCA embedding for diversity 
                                         dist_metric = "euclidean"
      )
      
    }
  )
  
  heatmap_func <- function(){
    quickPipeline(pathway_genes = v$pathway_genes, 
                  seurat_obj = master_seurat, 
                  k_final = v$opt_k, 
                  min_genes_on = input$min_genes_pathway, 
                  min_expr = input$min_expr_threshold
    )$p
  }
  
  heatmap_data_func <- function(){
    quickPipeline(pathway_genes = v$pathway_genes, 
                  seurat_obj = master_seurat, 
                  k_final = v$opt_k, 
                  min_genes_on = input$min_genes_pathway, 
                  min_expr = input$min_expr_threshold)$data_frame
  }
  
  output$heatmap <- renderPlot({
    
    input$ready
    
    isolate(
      heatmap_func()
    )
    
  })
  
  silhouette_func <- function(){
    v$silh_z
  }
  
  output$silh_z <- renderPlot({
    
    input$ready
    
    isolate(
      silhouette_func()
    )
    
  })
  
  output$d_silh_z <- downloadHandler(
    filename = function() {
      paste('silhouette_plot_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plot = silhouette_func(), device = "pdf")
    }
  )
  
  output$d_heatmap <- downloadHandler(
    filename = function() {
      paste('profile_heatmap_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plot = heatmap_func(), device = "pdf")
    }
  )
  
  output$d_heatmap_data <- downloadHandler(
    filename = function() {
      paste('profile_heatmap_data', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(file = file, x = heatmap_data_func())
    }
  )
  
  ecdf_diversity_func <- function(){
    ecdf_diversity(
      v$control_res)
  }
  
  output$ecdf_diversity <- renderPlot({
    
    input$ready
    
    isolate(
      ecdf_diversity_func()
    )
  })
  
  output$d_ecdf_diversity <- downloadHandler(
    filename = function() {
      paste('d_ecdf_diversity_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plot = ecdf_diversity_func(), device = "pdf")
    }
  )
  
  rank_diversity_func <- function(){
    rank_diversity(pathway_genes = v$pathway_genes,
                   make_plot = T,
                   k_final = v$opt_k,
                   min_expr = input$min_expr_threshold,
                   min_genes_on = input$min_genes_pathway,
                   manual_embedding = pca_proj,
                   seurat_obj=master_seurat
    )
  }
  
  output$rank_diversity <- renderPlot({
    
    input$ready
    
    isolate(
      rank_diversity_func()
    )
  })
  
  output$d_rank_diversity <- downloadHandler(
    filename = function() {
      paste('d_rank_diversity_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plot = rank_diversity_func(), device = "pdf")
    }
  )
  
  motif_heatmap_func <- function(){
    motif_heatmap(control_res=v$control_res,
                  pathway_genes=v$pathway_genes,
                  diverse_quantile=0.90,
                  type="motif"
    )
  }
  
  output$motif_heatmap <- renderPlot({
    
    input$ready
    
    isolate(
      motif_heatmap_func()
    )
  })
  
  output$d_motif_heatmap <- downloadHandler(
    filename = function() {
      paste('d_motif_heatmap_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)
      motif_heatmap_func()
      dev.off()
    }
  )
  
  motif_ct_heatmap_func <- function(){
    motif_ct_heatmap(control_res=v$control_res,
                     pathway_genes=v$pathway_genes,
                     diverse_quantile=0.90,
                     type="motif"
    )
  }
  
  output$motif_ct_heatmap <- renderPlot({
    
    input$ready
    
    isolate(
      motif_ct_heatmap_func()
    )
  })
  
  output$d_motif_ct_heatmap <- downloadHandler(
    filename = function() {
      paste('d_motif_ct_heatmap_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)
      motif_ct_heatmap_func()
      dev.off()
    }
  )
  
  glob_dendr_diversity_func <- function(){
    global_dendr(control_res = v$control_res, 
                 seurat_obj = master_seurat, 
                 use_pca = T, 
                 n_pcs = 1:20,
                 clust_method = 'ward.D2',
                 dist_metric ='cosine'
    )[[1]]
  }
  
  output$glob_dendr_diversity <- renderPlot({
    
    input$ready
    
    isolate(
      glob_dendr_diversity_func()
    )
  })
  
  output$d_glob_dendr_diversity <- downloadHandler(
    filename = function() {
      paste('d_glob_dendr_diversity_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      pdf(file, width = 4, height = 20)
      glob_dendr_diversity_func()
      dev.off()
    }
  )
  
  global_umap_func <- function(){
    global_umap(control_res = v$control_res, 
                seurat_obj = master_seurat, 
                use_pca = T, 
                n_pcs = 1:20,
                clust_method = 'ward.D2',
                dist_metric ='cosine'
    )
  }
  
  output$global_umap <- renderPlot({
    
    input$ready
    
    isolate(
      global_umap_func()
    )
  })
  
  output$d_global_umap <- downloadHandler(
    filename = function() {
      paste('d_global_umap_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)
      global_umap_func()
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)