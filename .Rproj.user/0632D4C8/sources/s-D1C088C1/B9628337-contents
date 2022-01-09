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
    # title = "A cellular atlas of DNA methylation dysregulation during mouse early organogenesis",
    title = "foo",
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
      sidebarPanel(width=3, selectizeInput(inputId = "gene", label = "Select gene", choices=NULL, selected="T")
      ),
      mainPanel(
        girafeOutput("plot_gene_expr_pseudoubulk")
      )
    ),

        
    tabPanel(
      title = "Differential expression (pseudobulk)", id = "diff_expr_pseudoubulk",
      sidebarPanel(width=3,
                   selectInput(inputId = "classA", label = "Select class A", choices = classes, selected = "WT"),
                   selectInput(inputId = "classB", label = "Select class B", choices = classes, selected = "Dnmt1_KO"),
                   plotOutput("differential_motif", height="100px", width="300px"),
                   # sliderInput("differential_range", label = "Range of differential values", min=-25, max=25, step=0.5, value = c(-25,25)),
                   # sliderInput("highlight_top_n_genes", label = "Highlight top N genes", min=0, max=100, step=1, value = 0)
      ),
      ),
      mainPanel(
        girafeOutput("plot_diff_expr_pseudoubulk")
      )
    ),
    
    tabPanel(
      title = "Celltype proportions", id = "celltype_proportions",
      sidebarPanel(width=3, 
        selectizeInput(inputId = "class", label = "Select class", choices=classes),
        checkboxInput("split_samples", "Split samples")
      ),
      mainPanel(
        girafeOutput("plot_celltype_proportions")
      )
    )
  )
)
