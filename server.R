library(R.utils)
library(shiny)
library(shinyFiles)
library(ggplot2)
library(HDF5Array)
library(data.table)
library(purrr)
library(cowplot)
library(ggrepel)
library(GGally)
library(ggiraph)
require(patchwork)
require(ggpubr) # to remove this dependency?
# library(DT)
# library(plotly)

setwd("/Users/argelagr/shiny_dnmt_tet")

basedir <- "/Users/argelagr/data/shiny_dnmt_tet"

################
## Load utils ##
################

source("utils.R")

# Updated ggnet2 function
source("/Users/argelagr/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/pseudobulk/ggnet2.R")

###############
## Load data ##
###############

# source("load_data.R")

###############
## Shiny app ##
###############

server <- function(input, output, session) {
  
  
  #######################
  ## Selectize speedup ##
  #######################
  
  updateSelectizeInput(session = session, inputId = 'gene', choices = genes, server = TRUE, selected = "T") 
  
  ##########
  ## UMAP ##
  ##########
  
  plot_UMAP <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$class <- "WT"
    # input$colourby <- "celltype"
    ## END TEST ##
    
    ## Fetch data ##
    
    # Select UMAP
    umap.dt <- umap_list[[input$class]]
    to.plot <- umap.dt %>% merge(sample_metadata[,c("cell","celltype","dataset","sample","stage")], by=c("cell"))
    
    # define color variable
    if (input$colourby == "gene_expression") {
      tmp <- data.table(
        cell = colnames(link_rna_expr),
        color = as.numeric(link_rna_expr[input$gene_umap_rna,])
      )
      to.plot <- to.plot %>% merge(tmp,by="cell") %>% setorder(color)
    } else {
      to.plot$color <- to.plot[[input$colourby]]
    }

    ## Plot PAGA ##
    
    celltypes <- sapply(paga$val,"[[","vertex.names")
    alphas <- rep(1.0,length(celltypes)); names(alphas) <- celltypes
    sizes <- rep(8,length(celltypes)); names(sizes) <- celltypes
    
    p.paga <- ggnet2_interactive(
      net = paga,
      mode = c("x", "y"),
      color = celltype_colours[celltypes],
      node.alpha = alphas,
      node.size = sizes,    
      edge.size = 0.15,
      edge.color = "grey",
      label = TRUE,
      label.size = 3
    )
    
    ## Plot UMAP ##
    
    # Subset UMAP for faster visualisation
    # to.plot <- to.plot[sample(1:.N, size = 10000)]
    
    p.umap <- ggplot(to.plot, aes(x = UMAP1, y = UMAP2, color = color)) +
      # geom_point(size = 1, alpha = 0.9) +
      geom_point_interactive(aes(tooltip = celltype, data_id = celltype), size=1, alpha=0.9) +
      coord_fixed(ratio = 0.8) +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "none"
      )

    # Modify legends    
    if (input$colourby%in%c("celltype","sample","dataset")) {
      p.umap <- p.umap + 
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
    } else {
      p.umap <- p.umap +
        theme(
          legend.title = element_blank()
        )
    }
    
    # Define palette
    palette <- switch(input$colourby, 
      "celltype" = celltype_palette, 
      "sample" = sample_palette, 
      "dataset" = dataset_palette, 
      "gene_expression" = rna_palette
    )
    p.umap <- p.umap + palette
      
    ## Barplot/Violin plot with statistics ##
    
    if (input$colourby%in%c("celltype","dataset","sample")) {
      
      if (input$colourby=="celltype") {
        to.plot2 <- to.plot[,.N,by=c("sample","celltype")]
        x_axis <- "sample"
      } else {
        to.plot2 <- to.plot[,.N,by=c(input$colourby,"celltype")]
        x_axis <- input$colourby
      }
      
      p2 <- ggplot(to.plot2, aes_string(x = x_axis, y = "N", fill = "celltype")) +
        geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat = "identity") +
        labs(y = "Number of cells") +
        celltype_palette_fill +
        theme_classic() +
        theme(
          legend.position = "none",
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black", size=rel(0.8)),
          axis.title.x = element_blank()
        )
      
    } else if (input$colourby=="gene_expression") {
      
      clust.sizes <- table(to.plot$celltype)
      
      p2 <- ggplot(to.plot, aes(x = celltype, y = color, fill = celltype)) +
        geom_violin(scale = "width", alpha=0.8) +
        # geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        geom_boxplot_interactive(aes(tooltip=celltype, data_id=celltype), width=0.5, outlier.shape=NA, alpha=0.8) +
        labs(y = "Log2 normalised counts") +
        annotate(
          geom = "text",
          # x = factor(names(clust.sizes)),
          x = names(clust.sizes),
          y = rep_len(c(max(to.plot$color)*1.1, max(to.plot$color) * 1.2), length.out = length(clust.sizes)),
          label = as.vector(clust.sizes),
          size = 3
        ) +
        celltype_palette_fill +
        theme_classic() +
        theme(
          axis.title = element_text(size = 11, color="black"),
          axis.text.y = element_text(size = 9, color="black"),
          axis.text.x = element_text(size = 12, color="black", angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }

    layout <- "
    ABB
    ABB
    ABB
    ABB
    CCC
    "
    girafe(
      # code = print((p.paga+p.umap)/p2 + plot_layout(widths = c(1,10), heights=c(4,1))),
      code = print(p.paga+p.umap+p2 + plot_layout(design = layout)),
      width_svg = 13, height_svg = 9,
      options = list( 
        opts_sizing(rescale = FALSE),
        # opts_selection(type = "single", css = "cursor:pointer;fill:magenta;color:magenta"),
        opts_selection(type = "single", css = ""),
        # opts_hover_inv(css = "opacity:0.45;"),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
  
  })
  
  output$umap = renderGirafe({
    # shiny::validate(need(input$gene_umap_rna%in%genes, "" ))
    plot_UMAP()
  })
 
  

  ##################################
  ## Gene expression (pseudobulk) ##
  ##################################
  
  ############################################
  ## Differential expression: volcano plots ##
  ############################################
  
  
  # p <- ggplot(to.plot, aes(x=celltype, y=class, fill=diff)) +
  #   geom_tile(color="black") +
  #   # scale_fill_gradientn(colours = terrain.colors(10), na.value = 'gray70') +
  #   # scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
  #   theme_classic() +
  #   guides(x = guide_axis(angle = 90)) +
  #   theme(
  #     axis.text = element_text(color="black", size=rel(0.7)),
  #     axis.title = element_blank(),
  #     strip.background = element_blank(),
  #     axis.ticks = element_blank(),
  #     axis.line = element_blank(),
  #     legend.title = element_blank()
  #   )
  
  ##########################
  ## Celltype proportions ##
  ##########################
  
  plot_celltype_proportions <- reactive({
    
    # Boxplots if N>1 and group by class, unless barplot option selected
    # Barplots if N==1 and group by class
    # Barplots if N==1 and group by sample
    
    ## START TEST ##
    input <- list()
    input$class <- "WT"
    input$colourby <- "celltype"
    input$split_samples <- FALSE
    ## END TEST ##
    
    to.plot <- sample_metadata %>% 
      .[class==input$class] %>%
      .[,N:=.N,by="sample"] %>%
      .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype")]

    # Define colours and cell type order
    tmp <- celltype_colours[names(celltype_colours) %in% unique(to.plot$celltype)]
    to.plot[,celltype:=factor(celltype, levels=rev(names(celltype_colours)))]

    samples.to.plot <- unique(to.plot$sample)
    
    # Plot
    
    if (input$split_samples) {
      p <- ggplot(to.plot, aes(x=celltype, y=N)) +
        geom_bar(aes(fill=celltype), stat="identity", color="black") +
        scale_fill_manual(values=celltype_colours) +
        facet_wrap(~sample, scales="fixed") +
        coord_flip() +
        labs(y="Number of cells") +
        theme_bw() +
        theme(
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(color="black", size=rel(0.9)),
          axis.title.x = element_text(color="black", size=rel(0.9)),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=rel(1), color="black"),
          axis.text.x = element_text(size=rel(1), color="black")
        )
    } else {
      p <- ggplot(to.plot, aes(x=celltype, y=N)) +
        geom_boxplot(aes(fill=celltype), color="black", outlier.shape=NA) +
        scale_fill_manual(values=celltype_colours) +
        coord_flip() +
        labs(y="Number of cells") +
        theme_bw() +
        theme(
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(color="black", size=rel(0.9)),
          axis.title.x = element_text(color="black", size=rel(0.9)),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=rel(1), color="black"),
          axis.text.x = element_text(size=rel(1), color="black")
        )
    }
    
    girafe(
      code = print(p),
      width_svg = 13, height_svg = 9,
      options = list( 
        opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
    
  })
  
  output$plot_celltype_proportions = renderGirafe({
    # shiny::validate(need(input$gene_umap_rna%in%genes, "" ))
    plot_celltype_proportions()
  })


}