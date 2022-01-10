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
  
  updateSelectizeInput(session = session, inputId = 'gene_umap_rna', choices = genes, server = TRUE, selected = "T") 
  updateSelectizeInput(session = session, inputId = 'gene_pseudobulk', choices = genes, server = TRUE, selected = "T") 
  
  ##########
  ## UMAP ##
  ##########
  
  plot_UMAP <- reactive({
    
    # TO-DO:
    # - subset number of cells
    # - add legend for sample
    # - add option to select samples and data sets
    # - add option to select celltypes
    # - increase size of UMAP or decrease size of cells when colouring by gene expr
    
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
  
  # TO-DO:
  # - fix order of classes 
  # - fix colors
  # - fix facet text
  # - add option to select celltypes
  # - add option to add number of cells for each barplot
  
  plot_gene_expr_pseudobulk <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$gene_pseudobulk <- "T"
    # input$dataset_gene_expr_pseudoubulk <- c("KO", "CRISPR")
    ## END TEST ##
    
    to.plot <- data.table(
      sample = colnames(sce.pseudobulk),
      expr = logcounts(sce.pseudobulk)[input$gene_pseudobulk,],
      class = sce.pseudobulk$class,
      celltype = sce.pseudobulk$celltype,
      dataset = sce.pseudobulk$dataset
    ) %>% .[dataset%in%input$dataset_gene_expr_pseudoubulk]
    
    p_list <- list()
    for (j in unique(to.plot$dataset)) {
      
      p_list[[j]] <- ggplot(to.plot[dataset==j], aes(x=class, y=expr, fill=class)) +
        geom_bar(stat="identity", color="black", width=0.75) +
        facet_wrap(~celltype, scales="fixed") +
        scale_fill_brewer(palette="Dark2") +
        # scale_fill_manual(values=opts$classes.colors) +
        theme_classic() +
        labs(x="",y="RNA expression", title=j) +
        guides(x = guide_axis(angle = 90)) +
        theme(
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size=rel(0.85)),
          axis.text.x = element_text(colour="black",size=rel(0.9)),
          axis.text.y = element_text(colour="black",size=rel(0.9)),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(colour="black",size=rel(1.0)),
          legend.position = "none"
          # legend.title = element_blank(),
          # legend.text = element_text(size=rel(0.85))
        )
    }
  return(cowplot::plot_grid(plotlist=p_list, ncol = 2))
  })
    
  output$plot_gene_expr_pseudoubulk = renderPlot({
    shiny::validate(need(input$gene_pseudobulk%in%genes, "" ))
    plot_gene_expr_pseudobulk()
  })
    
    
  ############################################
  ## Differential expression: volcano plots ##
  ############################################
  
  
  plot_diff_expr_volcano <- reactive({
  
    ## START TEST ##
    # input <- list()
    # input$classA_diff_expr_volcano <- "Dnmt1_KO"
    # input$classB_diff_expr_volcano <- "WT"
    # input$celltypes_diff_expr_volcano <- c("NMP","Spinal_cord")
    ## END TEST ##
    
    to.plot <- input$celltypes_diff_expr_volcano %>% map(function(i) {
      fread(file.path(basedir,sprintf("differential/Dnmt1_KO/%s_%s_vs_%s.txt.gz",i,input$classB_diff_expr_volcano,input$classA_diff_expr_volcano))) %>%
        .[,celltype:=i]
    }) %>% rbindlist %>% .[!is.na(logFC)]
    
    to.plot %>% .[, sig := (padj_fdr<=0.01 & abs(logFC)>=1)]    

    # top_genes <- 15
    
    # negative_hits <- to.plot[sig==TRUE & logFC<0,gene]
    # positive_hits <- to.plot[sig==TRUE & logFC>0,gene]
    # all <- nrow(to.plot)
  
    xlim <- max(abs(to.plot$logFC), na.rm=T)
    ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)
  
    p <- ggplot(to.plot, aes(x=logFC, y=-log10(padj_fdr+1e-100))) +
      labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)"))) +
      # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
      geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
      ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
      scale_color_manual(values=c("black","red")) +
      scale_size_manual(values=c(0.75,1.25)) +
      # scale_x_continuous(limits=c(-xlim-1.5,xlim+1.5)) +
      # scale_y_continuous(limits=c(0,ylim+15)) +
      # annotate("text", x=0, y=ylim+14, size=4, label=sprintf("(%d)", all)) +
      # annotate("text", x=-xlim-0.5, y=ylim+14, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      # annotate("text", x=xlim+0.5, y=ylim+14, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      # annotate("text", x=-xlim, y=0, size=3, label=sprintf("Up in %s (N=%s)",groupA,to.plot[["groupA_N"]][[1]])) +
      # annotate("text", x=xlim, y=0, size=3, label=sprintf("Up in %s (N=%s)",groupB,to.plot[["groupB_N"]][[1]])) +
      # ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=top_genes), aes(x=logFC, y=-log10(padj_fdr+1e-100), label=gene), size=3, max.overlaps = Inf) +
      facet_wrap(~celltype) +
      theme_classic() +
      theme(
        strip.background = element_blank(),
        axis.text = element_text(size=rel(0.75), color='black'),
        axis.title = element_text(size=rel(1.0), color='black'),
        legend.position="none"
      )
    return(p)
  })
  
  output$plot_diff_volcano = renderPlot({
    # shiny::validate(need(input$gene_umap_rna%in%genes, "" ))
    plot_diff_expr_volcano()
  })
  

  
  plot_diff_heatmap <- reactive({
    
    # TO-DO:
    # - split by dataset
    # - click on a gene allows you to see the barplots/boxplots
    
    
    ## START TEST ##
    # input <- list()
    # input$genes_diff_heatmap <- c("T","Sox2","Foxa2","Hoxd9")
    # input$celltypes_diff_heatmap <- celltypes
    # input$classes_diff_heatmap <- classes[classes!="WT"]
    ## END TEST ##
    
    to.plot <- expand.grid(input$celltypes_diff_heatmap, input$classes_diff_heatmap, input$genes_diff_heatmap) %>% 
      as.data.table %>% setnames(c("celltype","class","gene")) %>%
      merge(diff_pseudobulk.dt[gene%in%input$genes_diff_heatmap & celltype%in%input$celltypes_diff_heatmap & class%in%input$classes_diff_heatmap] , by=c("celltype","class","gene"), all.x=T) %>%
      .[,class:=factor(class,levels=input$classes_diff_heatmap)]
    
    # filter entries with lots of NAs
    to.plot <- to.plot[,foo:=mean(is.na(diff)),by=c("celltype","gene")] %>% .[foo<1] %>% .[,foo:=NULL]
    
    p <- ggplot(to.plot, aes(x=celltype, y=class, fill=diff)) +
      geom_tile(color="black") +
      facet_wrap(~gene) +
      # scale_fill_gradientn(colours = terrain.colors(10), na.value = 'gray70') +
      # scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        axis.text = element_text(color="black", size=rel(0.75)),
        axis.title = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank()
      )
    
    return(p)
    
  })
  
  output$plot_diff_heatmap = renderPlot({
    # shiny::validate(need(input$gene_umap_rna%in%genes, "" ))
    plot_diff_heatmap()
  })
  
  
  ##########################
  ## Celltype proportions ##
  ##########################
  
  plot_celltype_proportions <- reactive({
    
    print(input)
    
    ## START TEST ##
    # input <- list()
    # input$class <- "WT"
    # input$dataset <- c("KO", "CRISPR")
    # input$split_samples <- TRUE
    ## END TEST ##
    
    to.plot <- sample_metadata %>% 
      .[class==input$class_celltype_proportions & dataset%in%input$dataset_celltype_proportions] 
    
    # remove ExE cells
    if (input$remove_extraembryonic) {
      to.plot <- to.plot[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
    }
    
    # calculate celltype proportions
    to.plot <- to.plot %>%
      .[,N:=.N,by="sample"] %>%
      .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype","dataset")]

    # Define colours and cell type order
    tmp <- celltype_colours[names(celltype_colours) %in% unique(to.plot$celltype)]
    to.plot[,celltype:=factor(celltype, levels=rev(names(celltype_colours)))]
    to.plot[,sample:=sprintf("(%s) %s",dataset,sample)]
    
    # Plot
    if (input$split_samples) {
      p <- ggplot(to.plot, aes(x=celltype, y=celltype_proportion)) +
        # geom_bar(aes(fill=celltype), stat="identity", color="black") +
        geom_bar_interactive(aes(tooltip=celltype, data_id=celltype, fill=celltype), stat = "identity", color="black") +
        scale_fill_manual(values=celltype_colours) +
        facet_wrap(~sample, scales="free_x") +
        coord_flip() +
        labs(y="Fraction of cells") +
        theme_bw() +
        # guides(fill=guide_legend(ncol=1))
        guides(fill=guide_legend(ncol=1)) +
        theme(
          legend.position = "right",
          legend.key.size = unit(0.50, "cm"),
          # strip.background = element_blank(),
          # strip.text = element_text(color="black", size=rel(0.9)),
          axis.title.x = element_text(color="black", size=rel(0.75)),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=rel(0.5), color="black"),
          axis.text.x = element_text(size=rel(1), color="black")
        )
    } else {
      # stop()
      p <- ggplot(to.plot, aes(x=celltype, y=celltype_proportion, fill=celltype)) +
        geom_point_interactive(aes(tooltip=celltype, data_id=celltype), color="black", shape=21) +
        geom_boxplot_interactive(aes(tooltip=celltype, data_id=celltype), color="black", outlier.shape=NA, alpha=0.8) +
        scale_fill_manual(values=celltype_colours) +
        facet_wrap(~dataset, scales="fixed") +
        coord_flip() +
        labs(y="Fraction of cells") +
        theme_bw() +
        theme(
          legend.position = "none",
          strip.text = element_text(color="black", size=rel(1.25)),
          axis.title.x = element_text(color="black", size=rel(1.0)),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=rel(1.25), color="black"),
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