library(data.table)
library(purrr)
library(shiny)
library(shinyFiles)
library(ggplot2)
library(DT)
library(ggiraph)
library(rintrojs)
library(shinythemes)

#####################
## Define settings ##
#####################

basedir <- "/Users/argelagr/data/shiny_dnmt_tet"

source("utils.R")

big_plot_width = "900px"
big_plot_height = "500px"

narrower_plot_width = "650px"

half_plot_width = "450px"
narrower_half_plot_width = "350px"
half_plot_height = "260px"

###############
## Load data ##
###############

# celltypes <- fread(paste0(basedir,"/celltypes.txt"), header=F)[[1]]
genes <- fread(paste0(basedir,"/genes.txt"), header=F)[[1]]

##################
## Shiny app UI ##
##################

ui <- shinyUI(fluidPage(
  # sidebarLayout(
  
  # mainPanel(
  #   id = "main",
  #   width = 10,
  #   titlePanel(
  #     "A cellular atlas of DNA methylation dysregulation during mouse early organogenesis"
  #   ),
  navbarPage(
  # tabsetPanel(
    # id = "tabs",
    title = "A cellular atlas of DNA methylation dysregulation during mouse early organogenesis",
    # title = "foo",
    theme = shinytheme("spacelab"),
    
    tabPanel(
      title = "UMAP", id = "umap",
      sidebarPanel(width=3,
        selectInput(inputId = "class", label = "Class", choices = classes, selected = "WT"),
        selectInput(inputId = "colourby", label = "Plot colour", choices = c("Cell type"="celltype", "Data set"="dataset", "Sample"="sample", "Gene expression"="gene_expression"), selected = "celltype"),
        conditionalPanel(
          condition = "input.colourby == 'gene_expression'",
          selectizeInput("gene_umap_rna", "Select gene to show RNA expression", choices = NULL, selected = "T")
        )
        # checkboxInput("numbers", "Annotate clusters in plot"),
      ),
      mainPanel(
        girafeOutput("umap", width = "900px", height = "800px"),
        # plotOutput("stage_contribution", width = big_plot_width)
      )
    ),
    
    tabPanel(
      title = "Gene expression (pseudobulk)", id = "gene_expr_pseudoubulk",
      sidebarPanel(width=3, 
        selectizeInput(inputId = "gene_pseudobulk", label = "Select gene", choices=NULL, selected="T"),
        selectInput("classes_gene_expr_pseudobulk", "Classes", choices = classes, selected = classes, multiple = TRUE),
        selectInput("celltypes_gene_expr_pseudobulk", "Celltypes", choices = celltypes, selected = celltypes, multiple = TRUE),
        checkboxGroupInput("dataset_gene_expr_pseudoubulk", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR"))
      ),
      mainPanel(
        plotOutput("plot_gene_expr_pseudoubulk",  width = "1050px", height = "800px")
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
                   selectInput("celltypes_diff_heatmap", "Celltypes", choices = celltypes, selected = celltypes, multiple = TRUE)
      ),
      mainPanel(
        HTML("A positive sign indicates that the gene is more expressed in the WT"),
        plotOutput("plot_diff_heatmap")
      )
    ),
    
    tabPanel(
      title = "Celltype proportions", id = "celltype_proportions",
      sidebarPanel(width=3, 
        selectizeInput(inputId = "class_celltype_proportions", label = "Select class", choices=classes, selected="WT"),
        checkboxInput("split_samples", "Split samples", value = T),
        checkboxInput("remove_extraembryonic", "Remove Extraembryonic cell types", value = F),
        checkboxGroupInput("dataset_celltype_proportions", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR"))
      ),
      mainPanel(
        girafeOutput("plot_celltype_proportions")
      )
    ),
    
    tabPanel(
      title = "Celltype proportions comparisons", id = "celltype_proportions_comparisons",
      sidebarPanel(width=3, 
                   selectizeInput(inputId = "class_celltype_comparisons", label = "Select class", choices=classes, selected="DNMT1_KO"),
                   checkboxInput("split_samples_celltype_comparisons", "Split samples", value = T),
                   checkboxInput("remove_extraembryonic_celltype_comparisons", "Remove Extraembryonic cell types", value = F),
                   checkboxGroupInput("dataset_celltype_comparisons", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR")),
                   selectInput("output_type_celltype_comparisons", label = "Output type",  choices = list("Box plots" = "box_plots", "Polar plots" = "polar_plots"), selected = "Boxplots")
      ),
      mainPanel(
        # girafeOutput("plot_celltype_proportions_comparisons")
        plotOutput("plot_celltype_proportions_comparisons")
      )
    ),
    
    
    
    tags$style(HTML(".navbar-header { width:100% } .navbar-brand { width: 100%; text-align: left; font-size: 150%; }"))
  )
))
