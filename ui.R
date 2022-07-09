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
    
    #############
    ## Mapping ##
    #############
    
    tabPanel(
      title = "Mapping", id = "mapping",
      sidebarPanel(width=3, 
                   selectizeInput("mapping_sample", "Select sample", choices = samples, selected = "WT_1"),
                   checkboxInput("mapping_subset_cells", "Subset atlas cells", value=TRUE)
                   # selectInput("mapping_sample", "Select sample", choices = samples, selected = samples, multiple = TRUE),
      ),
      mainPanel(
        girafeOutput("mapping_plot", width = "1000px", height = "600px")
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
      title = "Diff. expr (volcano plots)", id = "diff_volcano",
      sidebarPanel(width=3,
                   selectInput(inputId = "diff_volcano_class", label = "Select class", choices = classes[classes!="WT"], selected = "Dnmt1_KO"),
                   selectInput(inputId = "diff_resolution", label = "Select resolution", choices = c("Cells","Pseudobulk"), selected = "Pseudobulk"),
                   selectInput("diff_volcano_celltype", "Celltypes", choices = celltypes, selected = "NMP", multiple = TRUE),
                   selectizeInput("diff_gene", "Select feature", choices = NULL),
                   sliderInput("diff_min_log_pval", label = "Minimum log pvalue", min=0, max=100, step=1, value = 25),
                   sliderInput("diff_range", label = "Range of differential values", min=-8, max=8, step=0.5, value = c(-8,8)),
      ),
      mainPanel(
        # HTML("A positive sign indicates that the gene is more expressed in the WT"),
        plotOutput("plot_diff_volcano")
      )
    ),
    
    tabPanel(
      title = "Differential analysis", id = "diff",
      sidebarPanel(width=3,
                   selectInput(inputId = "diff_celltypeA", label = "Select celltype A", choices = celltypes, selected = "Gut"),
                   selectInput(inputId = "diff_celltypeB", label = "Select celltype B", choices = celltypes, selected = "Neural_crest"),
                   
                   selectInput(inputId = "diff_modality", label = "Select data modality", choices = c("RNA","ATAC"), selected = "RNA"),
                   conditionalPanel(
                     condition = "input.diff_modality == 'ATAC'",
                     selectInput("diff_atac_chr", "Select chromosome", choices = chr_mm10, selected = "chr1")
                   ),

                   # sliderInput("highlight_top_n_genes", label = "Highlight top N genes", min=0, max=100, step=1, value = 0)
      ),
      mainPanel(
        girafeOutput("diff_plot", width = "900px", height = "600px"),
        HTML("For visualisation efficiency, there is a maximum of 10,000 features to be displayed")
      )
    ),
    
    
    tabPanel(
      title = "Diff. expr (heatmap)", id = "diff_heatmap_pseudoubulk",
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
