# for testing
# shiny::loadSupport()

#############################
## Define global variables ##
#############################

# celltypes <- fread(paste0(data_folder,"/celltypes.txt"), header=F)[[1]]
genes <- fread(paste0(data_folder,"/rna_expression/genes.txt"), header=F)[[1]]
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
    
    #############
    ## Mapping ##
    #############
    
    tabPanel(
      title = "Mapping", id = "mapping",
      sidebarPanel(width=3, 
       selectizeInput("mapping_sample", "Select sample", choices = samples, selected = "WT_1"),
       checkboxInput("mapping_subset_cells", "Subset atlas cells", value=TRUE)
      ),
      mainPanel(
        girafeOutput("mapping_plot", width = "1000px", height = "600px")
      )
    ),
    
    ##########
    ## UMAP ##
    ##########
    
    tabPanel(
      title = "UMAP", id = "umap",
      sidebarPanel(width=3,
        selectInput(inputId = "class", label = "Class", choices = classes, selected = "WT"),
        selectInput(inputId = "colourby", label = "Plot colour", choices = c("Cell type"="celltype", "Data set"="dataset", "Sample"="sample", "Gene expression"="gene_expression"), selected = "celltype"),
        checkboxInput("subset_cells_umap", "Subset number of cells for the UMAP", value = TRUE),
        conditionalPanel(
          condition = "input.colourby == 'gene_expression'",
          selectizeInput("gene_umap", "Select gene to show RNA expression", choices = NULL, selected = "T")
        )
        # checkboxInput("numbers", "Annotate clusters in plot"),
      ),
      mainPanel(
        girafeOutput("umap", width = "1000px", height = "800px"),
        # plotOutput("stage_contribution", width = big_plot_width)
      )
    ),
    
    #####################
    ## Gene expression ##
    #####################
    
    tabPanel(
      title = "Gene expression", id = "gene_expr",
      sidebarPanel(width=3, 
        selectizeInput(inputId = "gene_expr_gene", label = "Select gene", choices=NULL, selected="T"),
        selectInput("gene_expr_classes", "Classes", choices = classes, selected = classes, multiple = TRUE),
        selectInput("gene_expr_celltypes", "Celltypes", choices = celltypes_pseudobulk, selected = celltypes_pseudobulk, multiple = TRUE),
        checkboxGroupInput("gene_expr_dataset", label = "Data set",  choices = list("KO" = "KO", "CRISPR (Grosswendt2020)" = "CRISPR"), selected = c("KO","CRISPR")),
        selectInput("gene_expr_resolution", "Resolution", choices = c("Cells","Pseudobulk"), selected = "Pseudobulk"),
        checkboxInput("gene_expr_add_number_observations", "Show number of pseudobulk replicates per condition", value = TRUE)
      ),
      mainPanel(
        plotOutput("plot_gene_expr",  width = "1000px", height = "700px")
        # girafeOutput("plot_gene_expr")
      )
    ),

    #############################################
    ## Differential expression (volcano plots) ##
    #############################################
    
    tabPanel(
      title = "Diff. expr (volcano plots)", id = "diff_volcano",
      sidebarPanel(width=3,
       selectInput(inputId = "diff_class", label = "Select class", choices = classes[classes!="WT"], selected = "Dnmt1_KO"),
       selectInput(inputId = "diff_resolution", label = "Select resolution", choices = c("Cells"), selected = "Cells"),
       selectInput("diff_celltype", "Celltype", choices = celltypes, selected = "NMP"),
       selectizeInput("diff_gene", "Select gene", choices = NULL, selected = "Rhox5"),
       sliderInput("diff_min_log_pval", label = "Minimum log pvalue", min=0, max=100, step=1, value = 5),
       sliderInput("diff_range", label = "Range of differential values", min=-5, max=5, step=0.5, value = c(-5,5)),
      ),
      mainPanel(
        girafeOutput("diff_plot", width = "900px", height = "600px")
        # HTML("For visualisation efficiency, there is a maximum of 10,000 features to be displayed")
      )
    ),
    
    ########################################
    ## Differential expression (heatmaps) ##
    ########################################
    
    tabPanel(
      title = "Diff. expr (heatmap)", id = "diff_heatmap",
      sidebarPanel(width=3,
       selectInput("diff_heatmap_genes", "Genes", choices = NULL, selected = c("Hoxb9","Hoxc9","Hoxd9","Utf1","Slc7a3","Pim2","Xlr3a","Rhox9","Apoe"), multiple = TRUE),
       selectInput("diff_heatmap_classes", "Classes", choices = classes[classes!="WT"], selected = classes[classes!="WT"], multiple = TRUE),
       selectInput("diff_heatmap_celltypes", "Celltypes", choices = celltypes_pseudobulk, selected = celltypes_pseudobulk, multiple = TRUE),
       selectInput("diff_heatmap_split", "Facet by", choices = c("Celltype"="celltype", "Gene"="gene"), selected = "Gene"),
      ),
      mainPanel(
        HTML("A positive sign indicates that the gene is more expressed in the KO"),
        girafeOutput("plot_diff_heatmap")
      )
    ),
    
    ###########################
    ## Cell type proportions ##
    ###########################
    
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
    
    #########################
    ## Repetitive elements ##
    #########################
    
    tabPanel(
      title = "Repetitive elements", id = "repetitive",
      sidebarPanel(width=3,
        selectInput("repetitive_elements", "Elements", choices = repeat_classes, selected = c("IAP","LINE_L1","LTR_ERVK"), multiple = TRUE),
        selectInput("repetitive_classes", "Classes", choices = classes[classes!="WT"], selected = classes[classes!="WT"], multiple = TRUE),
        selectInput("repetitive_celltypes", "Celltypes", choices = celltypes_pseudobulk[1:4], selected = celltypes_pseudobulk, multiple = TRUE)
      ),
      mainPanel(
        HTML("A positive sign indicates that the repeat element is more expressed in the KO"),
        girafeOutput("plot_repetitive")
      )
    ),
    
    tags$style(HTML(".navbar-header { width:100% } .navbar-brand { width: 100%; text-align: left; font-size: 150%; }"))
  )
))
