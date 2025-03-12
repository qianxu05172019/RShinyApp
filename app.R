library(shiny)
library(Seurat)
library(DT)
library(ggplot2)
library(dplyr)

# UI definition
ui <- fluidPage(
  titlePanel("Seurat Analysis App"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("rds_file", "Upload Seurat object (.rds file)",
                accept = ".rds"),
      
      # Quality control parameters
      conditionalPanel(
        condition = "output.data_loaded",
        h4("Quality Control"),
        numericInput("min_features", "Minimum features per cell:", 200),
        numericInput("max_features", "Maximum features per cell:", 2500),
        numericInput("max_mt", "Maximum mitochondrial percentage:", 5),
        actionButton("run_qc", "Run QC")
      ),
      
      # Normalization and scaling
      conditionalPanel(
        condition = "output.qc_done",
        h4("Normalization"),
        selectInput("norm_method", "Normalization method:",
                   choices = c("LogNormalize", "SCTransform")),
        actionButton("run_norm", "Run Normalization")
      ),
      
      # Dimension reduction
      conditionalPanel(
        condition = "output.norm_done",
        h4("Dimension Reduction"),
        numericInput("n_pcs", "Number of PCs:", 30),
        numericInput("n_neighbors", "Number of neighbors:", 30),
        numericInput("min_dist", "Minimum distance:", 0.3),
        actionButton("run_dim_red", "Run Dimension Reduction")
      ),
      
      # Clustering
      conditionalPanel(
        condition = "output.dim_red_done",
        h4("Clustering"),
        numericInput("resolution", "Resolution:", 0.8),
        actionButton("run_cluster", "Run Clustering")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("QC Plots",
                 plotOutput("qc_violin"),
                 plotOutput("qc_scatter")),
        tabPanel("Feature Plots",
                 selectInput("feature_to_plot", "Select feature:", ""),
                 plotOutput("feature_plot")),
        tabPanel("Dimension Reduction",
                 selectInput("dim_red_type", "Plot type:",
                           choices = c("UMAP", "tSNE", "PCA")),
                 plotOutput("dim_red_plot")),
        tabPanel("Clustering Results",
                 plotOutput("cluster_plot")),
        tabPanel("Differential Expression",
                 selectInput("cluster_de", "Select cluster:", ""),
                 DT::dataTableOutput("de_table"))
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Reactive values to store state
  vals <- reactiveValues(
    seurat_obj = NULL,
    data_loaded = FALSE,
    qc_done = FALSE,
    norm_done = FALSE,
    dim_red_done = FALSE,
    cluster_done = FALSE
  )
  
  # Load data
  observeEvent(input$rds_file, {
    req(input$rds_file)
    vals$seurat_obj <- readRDS(input$rds_file$datapath)
    vals$data_loaded <- TRUE
    
    # Update feature selection
    updateSelectInput(session, "feature_to_plot",
                     choices = rownames(vals$seurat_obj))
  })
  
  # Run QC
  observeEvent(input$run_qc, {
    req(vals$seurat_obj)
    
    # Calculate percent mitochondrial
    vals$seurat_obj[["percent.mt"]] <- PercentageFeatureSet(vals$seurat_obj, 
                                                           pattern = "^MT-")
    
    # Filter cells
    vals$seurat_obj <- subset(vals$seurat_obj,
                             subset = nFeature_RNA > input$min_features &
                               nFeature_RNA < input$max_features &
                               percent.mt < input$max_mt)
    
    vals$qc_done <- TRUE
  })
  
  # Run normalization
  observeEvent(input$run_norm, {
    req(vals$seurat_obj, vals$qc_done)
    
    if(input$norm_method == "LogNormalize") {
      vals$seurat_obj <- NormalizeData(vals$seurat_obj)
    } else {
      vals$seurat_obj <- SCTransform(vals$seurat_obj)
    }
    
    # Find variable features
    vals$seurat_obj <- FindVariableFeatures(vals$seurat_obj)
    vals$seurat_obj <- ScaleData(vals$seurat_obj)
    
    vals$norm_done <- TRUE
  })
  
  # Run dimension reduction
  observeEvent(input$run_dim_red, {
    req(vals$seurat_obj, vals$norm_done)
    
    # Run PCA
    vals$seurat_obj <- RunPCA(vals$seurat_obj)
    
    # Run UMAP
    vals$seurat_obj <- RunUMAP(vals$seurat_obj,
                              dims = 1:input$n_pcs,
                              n.neighbors = input$n_neighbors,
                              min.dist = input$min_dist)
    
    # Run tSNE
    vals$seurat_obj <- RunTSNE(vals$seurat_obj,
                              dims = 1:input$n_pcs)
    
    vals$dim_red_done <- TRUE
  })
  
  # Run clustering
  observeEvent(input$run_cluster, {
    req(vals$seurat_obj, vals$dim_red_done)
    
    vals$seurat_obj <- FindNeighbors(vals$seurat_obj)
    vals$seurat_obj <- FindClusters(vals$seurat_obj,
                                   resolution = input$resolution)
    
    # Update cluster selection for DE
    updateSelectInput(session, "cluster_de",
                     choices = levels(vals$seurat_obj$seurat_clusters))
    
    vals$cluster_done <- TRUE
  })
  
  # QC plots
  output$qc_violin <- renderPlot({
    req(vals$seurat_obj)
    VlnPlot(vals$seurat_obj,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3)
  })
  
  output$qc_scatter <- renderPlot({
    req(vals$seurat_obj)
    FeatureScatter(vals$seurat_obj,
                   feature1 = "nCount_RNA",
                   feature2 = "nFeature_RNA")
  })
  
  # Feature plot
  output$feature_plot <- renderPlot({
    req(vals$seurat_obj, input$feature_to_plot)
    FeaturePlot(vals$seurat_obj,
                features = input$feature_to_plot)
  })
  
  # Dimension reduction plot
  output$dim_red_plot <- renderPlot({
    req(vals$seurat_obj, vals$dim_red_done)
    
    if(input$dim_red_type == "UMAP") {
      DimPlot(vals$seurat_obj, reduction = "umap")
    } else if(input$dim_red_type == "tSNE") {
      DimPlot(vals$seurat_obj, reduction = "tsne")
    } else {
      DimPlot(vals$seurat_obj, reduction = "pca")
    }
  })
  
  # Clustering plot
  output$cluster_plot <- renderPlot({
    req(vals$seurat_obj, vals$cluster_done)
    DimPlot(vals$seurat_obj, reduction = "umap", group.by = "seurat_clusters")
  })
  
  # Differential expression table
  output$de_table <- DT::renderDataTable({
    req(vals$seurat_obj, vals$cluster_done, input$cluster_de)
    
    markers <- FindMarkers(vals$seurat_obj,
                          ident.1 = input$cluster_de)
    
    DT::datatable(markers)
  })
  
  # Update conditional panel states
  output$data_loaded <- reactive({
    return(vals$data_loaded)
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  output$qc_done <- reactive({
    return(vals$qc_done)
  })
  outputOptions(output, "qc_done", suspendWhenHidden = FALSE)
  
  output$norm_done <- reactive({
    return(vals$norm_done)
  })
  outputOptions(output, "norm_done", suspendWhenHidden = FALSE)
  
  output$dim_red_done <- reactive({
    return(vals$dim_red_done)
  })
  outputOptions(output, "dim_red_done", suspendWhenHidden = FALSE)
}

# Run the app
shinyApp(ui = ui, server = server) 