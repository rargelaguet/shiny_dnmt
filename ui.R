# for testing
# shiny::loadSupport()

#############################
## Define global variables ##
#############################

celltypes <- fread(paste0(data_folder,"/celltypes.txt"), header=F)[[1]]
genes <- fread(paste0(data_folder,"/genes.txt"), header=F)[[1]]
classes <- c("WT","Dnmt3a_KO","Dnmt3b_KO","Dnmt1_KO")
celltypes_pseudobulk <- celltypes

##############
## Shiny UI ##
##############

ui <- shinyUI(fluidPage(
  navbarPage(
    title = "A cellular atlas of DNA methylation dysregulation during mouse early organogenesis",
    theme = shinytheme("spacelab"),

    tabPanel(
      title = "Overview", id = "overview",
      br(),
      includeMarkdown("overview.md")
    ),
    
    tabPanel(
      title = "Dataset stats", id = "stats",
      sidebarPanel(width=3,
        selectInput(inputId = "stat_to_plot", label = "Statistic to plot", choices = c("Number of cells"="ncells", "Number of embryos"="nembryos", "Number of genes"="ngenes"), selected = "ncells"),
      ),
      mainPanel(
        # girafeOutput("umap", width = "1000px", height = "800px"),
        plotOutput("dataset_stats", width = "800px", height = "400px")
      )
    ),
    
    tabPanel(
      title = "UMAP", id = "umap",
      sidebarPanel(width=3,
        selectInput(inputId = "class", label = "Class", choices = classes, selected = "WT"),
        selectInput(inputId = "colourby", label = "Plot colour", choices = c("Cell type"="celltype", "Data set"="dataset", "Sample"="sample", "Gene expression"="gene_expression"), selected = "celltype"),
        checkboxInput("subset_cells_umap", "Subset number of cells for the UMAP", value = TRUE),
        conditionalPanel(
          condition = "input.colourby == 'gene_expression'",
          selectizeInput("gene_umap_rna", "Select gene to show RNA expression", choices = NULL, selected = "T")
        )
        # checkboxInput("numbers", "Annotate clusters in plot"),
      ),
      mainPanel(
        girafeOutput("umap", width = "1000px", height = "800px"),
        # plotOutput("stage_contribution", width = big_plot_width)
      )
    ),
    
    tabPanel(
      title = "Gene expression", id = "gene_expr_pseudoubulk",
      sidebarPanel(width=3, 
        selectizeInput(inputId = "gene_pseudobulk", label = "Select gene", choices=NULL, selected="T"),
        selectInput("classes_gene_expr_pseudobulk", "Classes", choices = classes, selected = classes, multiple = TRUE),
        selectInput("celltypes_gene_expr_pseudobulk", "Celltypes", choices = celltypes_pseudobulk, selected = celltypes_pseudobulk, multiple = TRUE),
        checkboxGroupInput("dataset_gene_expr_pseudoubulk", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR")),
        checkboxInput("add_number_observations_pseudobulk", "Show number of pseudobulk replicates per condition", value = TRUE),
      ),
      mainPanel(
        plotOutput("plot_gene_expr_pseudoubulk",  width = "800px", height = "700px")
        # girafeOutput("plot_gene_expr_pseudoubulk")
      )
    ),

    tabPanel(
      title = "Diff. expr (volcano plots)", id = "diff_expr_volcano",
      sidebarPanel(width=3,
                   selectInput(inputId = "classA_diff_expr_volcano", label = "Select class A", choices = classes[classes!="WT"], selected = "Dnmt1_KO"),
                   selectInput(inputId = "classB_diff_expr_volcano", label = "Select class B", choices = "WT", selected = "WT"),
                   selectInput("celltypes_diff_expr_volcano", "Celltypes", choices = celltypes, selected = "NMP", multiple = TRUE)
      ),
      mainPanel(
        # HTML("A positive sign indicates that the gene is more expressed in the WT"),
        plotOutput("plot_diff_volcano")
      )
    ),
    
    
    tabPanel(
      title = "Diff. expr (heatmap)", id = "diff_expr_heatmap_pseudoubulk",
      sidebarPanel(width=3,
                   selectInput("genes_diff_heatmap", "Genes", choices = genes, selected = c("T","Sox2","Foxa2","Hoxd9"), multiple = TRUE),
                   selectInput("classes_diff_heatmap", "Classes", choices = classes[classes!="WT"], selected = classes[classes!="WT"], multiple = TRUE),
                   selectInput("celltypes_diff_heatmap", "Celltypes", choices = celltypes_pseudobulk, selected = celltypes_pseudobulk, multiple = TRUE)
      ),
      mainPanel(
        HTML("A positive sign indicates that the gene is more expressed in the WT"),
        # plotOutput("plot_diff_heatmap")
        girafeOutput("plot_diff_heatmap")
      )
    ),
    
    tabPanel(
      title = "Celltype proportions", id = "celltype_proportions",
      sidebarPanel(width=3, 
        selectizeInput(inputId = "class_celltype_proportions", label = "Select class", choices=classes, selected="WT"),
        # checkboxInput("split_samples", "Split samples", value = T),
        checkboxInput("remove_extraembryonic", "Remove Extraembryonic cell types", value = F),
        checkboxGroupInput("dataset_celltype_proportions", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR")),
        selectizeInput("visualisation_type_celltype_proportions", label = "visualisation",  choices = c("Barplots per sample", "Boxplots per class"), selected = "Boxplots per class")
      ),
      mainPanel(
        girafeOutput("plot_celltype_proportions")
      )
    ),
    
    tabPanel(
      title = "Celltype proportions comparisons", id = "celltype_proportions_comparisons",
      sidebarPanel(width=3, 
                   selectizeInput(inputId = "class_celltype_comparisons", label = "Select class", choices=classes[classes!="WT"], selected="Dnmt1_KO"),
                   # selectInput("class_celltype_comparisons", "Select class", choices = classes[classes!="WT"], selected = classes[classes!="WT"], multiple = TRUE),
                   checkboxInput("remove_extraembryonic_celltype_comparisons", "Remove Extraembryonic cell types", value = F),
                   checkboxInput("split_by_dataset_celltype_comparisons", label = "Split by data set",  value=F),
                   checkboxGroupInput("dataset_celltype_comparisons", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR"))
      ),
      mainPanel(
        girafeOutput("plot_celltype_proportions_comparisons")
      )
    ),
    
    
    
    tags$style(HTML(".navbar-header { width:100% } .navbar-brand { width: 100%; text-align: left; font-size: 150%; }"))
  )
))
