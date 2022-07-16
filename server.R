
# for testing
# shiny::loadSupport()
# source("load_data.R")

##################
## Shiny server ##
##################

server <- function(input, output, session) {
  
  ###################
  ## Dataset stats ##
  ###################
  
  plot_dataset_stats <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$stat_to_plot <- "ncells"
    ## END TEST ##
    
    if (input$stat_to_plot=="ncells") {
      to.plot <- cell_metadata.dt %>% .[,.N,by=c("class","dataset")]
      lab <- "Number of cells"
      plot_type <- "barplot"
    } else if (input$stat_to_plot=="nembryos") {
      to.plot <- cell_metadata.dt %>% 
        .[,.(N=length(unique(sample))), by=c("class","dataset")] %>%
        .[,N:=factor(N, levels=1:max(N))]
      lab <- "Number of embryos"
      plot_type <- "barplot"
    } else if (input$stat_to_plot=="ngenes") {
      to.plot <- cell_metadata.dt[,c("sample","class","dataset","nFeature_RNA")] %>% 
        setnames("nFeature_RNA","value")
      lab <- "Number of genes"
      plot_type <- "boxplot"
    }
    
    # c("ngenes","mit_percentage")
    
    if (plot_type=="barplot") {
      p <- ggbarplot(to.plot, x="class", y="N", fill="dataset", position=position_dodge(width = 0.75)) +
        labs(x="", y=lab) +
        scale_fill_brewer(palette="Dark2") +
        theme(
          legend.position = "top",
          legend.title = element_blank(),
          axis.text.y = element_text(colour="black",size=rel(0.8)),
          axis.text.x = element_text(colour="black",size=rel(1.0)),
        )
    } else if (plot_type=="boxplot") {
      
      to.plot.jitter <- to.plot %>% .[sample.int(n=nrow(.), size=nrow(.)/10)]
      
      p <- ggplot(to.plot, aes_string(x="sample", y="value", fill="class")) +
        ggrastr::geom_jitter_rast(aes(color=class), alpha=0.5, width=0.15, size=0.05, data=to.plot.jitter) +
        geom_boxplot(outlier.shape=NA, coef=1, alpha=0.9) +
        # facet_wrap(~variable, scales="free_y", labeller = as_labeller(facet.labels), nrow=1) +
        scale_fill_manual(values=class_colors) +
        scale_color_manual(values=class_colors) +
        guides(x = guide_axis(angle = 90)) +
        labs(x="", y=lab) +
        theme_classic() +
        theme(
          legend.title = element_blank(),
          legend.position = "top",
          axis.text.y = element_text(colour="black",size=rel(0.9)),
          # axis.text.x = element_text(colour="black",size=rel(0.55), angle=20, hjust=1, vjust=1),
          axis.text.x = element_text(colour="black",size=rel(0.85)),
          axis.title.x = element_blank()
        )
    }
    
    
    
    return(p)
    
  })
  
  output$dataset_stats = renderPlot({
    plot_dataset_stats()
  })
  
  
  #############
  ## Mapping ##
  #############
  
  plot_mapping <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$mapping_sample <- "WT_1"
    # input$mapping_subset_cells <- TRUE
    ## END TEST ##
    
    # Fetch data    
    cell_metadata_sample.dt <- cell_metadata.dt[sample==input$mapping_sample] %>% 
      .[,celltype:=stringr::str_replace_all(celltype,c("/"=" ","-"=" ","_"=" "))]
    
    # Define dot size  
    size.values <- c(0.30, 0.30)
    names(size.values) <- c(input$mapping_sample, "Atlas")
    
    # Define dot alpha  
    alpha.values <- c(0.70, 0.70)
    names(alpha.values) <- c(input$mapping_sample, "Atlas")
    
    # Define dot colours  
    colour.values <- c("red", "lightgrey")
    names(colour.values) <- c(input$mapping_sample, "Atlas")
    
    # Subset atlas for faster plotting
    if (input$mapping_subset_cells) {
      umap_reference_subset.dt <- rbind(
        umap_reference.dt[cell%in%unique(cell_metadata_sample.dt$closest.cell)],
        umap_reference.dt[!cell%in%unique(cell_metadata_sample.dt$closest.cell),.SD[sample.int(n=.N, size=1e4)]]
      )
    } else {
      umap_reference_subset.dt <- rbind(
        umap_reference.dt[cell%in%unique(cell_metadata_sample.dt$closest.cell)],
        umap_reference.dt[!cell%in%unique(cell_metadata_sample.dt$closest.cell),.SD[sample.int(n=.N, size=5e4)]]
      )
    }
    
    to.plot <- umap_reference_subset.dt %>% 
      .[,index:=match(cell, cell_metadata_sample.dt$closest.cell)] %>% 
      .[,mapped:=as.factor(!is.na(index))] %>% 
      .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",input$mapping_sample))] %>%
      setorder(mapped) 
    
    # p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas", rasterise = TRUE, subset = TRUE, legend = FALSE)
    
    p.mapping <- ggplot(to.plot, aes(x=V1, y=V2)) +
      geom_point_interactive(aes(size=mapped, alpha=mapped, colour=mapped, tooltip=celltype, data_id=celltype)) +
      scale_size_manual(values = size.values) +
      scale_alpha_manual(values = alpha.values) +
      scale_colour_manual(values = colour.values) +
      guides(colour = guide_legend(override.aes = list(size=6))) +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "top",
        legend.title = element_blank()
      )
    
    # Plot celltype proportions
    celltype_proportions.dt <- cell_metadata_sample.dt %>%
      .[,N:=.N,by="sample"] %>%
      .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype")]
    
    # tmp <- celltype_colours[names(celltype_colours) %in% unique(celltype_proportions.dt$celltype)]
    tmp <- celltype_colours; names(tmp) <- stringr::str_replace_all(names(tmp),c("/"=" ","-"=" ","_"=" "))
    # stopifnot(unique(cell_metadata_sample.dt$celltype)%in%names(tmp))
    # unique(cell_metadata_sample.dt$celltype)[!unique(cell_metadata_sample.dt$celltype)%in%names(tmp)]
    celltype_proportions.dt[,celltype:=factor(celltype, levels=rev(names(tmp)))]
    
    p.celltype_proportions <- ggplot(celltype_proportions.dt, aes(x=celltype, y=N)) +
      geom_bar_interactive(aes(tooltip=celltype, data_id=celltype, fill=celltype), stat="identity", color="black") +
      scale_fill_manual(values=tmp, drop=FALSE) +
      scale_x_discrete(drop=FALSE) +
      coord_flip() +
      labs(y="Number of cells") +
      theme_bw() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(color="black", size=rel(0.9)),
        axis.title.x = element_text(color="black", size=rel(0.9)),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=rel(0.9), color="black"),
        axis.text.x = element_text(size=rel(1), color="black")
      )
    
    girafe(
      code = print(p.mapping+p.celltype_proportions + plot_layout(nrow=1)),
      width_svg = 10, height_svg = 6,
      options = list(
        # opts_zoom(min = 1, max = 3),
        # opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;stroke:magenta;color:magenta;fill-opacity:1")
      )
    ) %>% return(.)
    
    
  })
  
  output$mapping_plot = renderGirafe({
    shiny::validate(need(input$mapping_sample%in%samples,""))
    plot_mapping()
  })
  
  ##########
  ## UMAP ##
  ##########
  
  updateSelectizeInput(session = session, inputId = 'gene_umap', choices = genes, server = TRUE, selected = "T") 
  
  plot_UMAP <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$class <- "WT"
    # input$colourby <- "celltype"
    ## END TEST ##
    
    ## Fetch data ##
    
    # Load UMAP
    umap.dt <- fread(file.path(data_folder,sprintf("/dimensionality_reduction/%s/umap.txt.gz",input$class)))
    to.plot <- umap.dt %>% merge(cell_metadata.dt[,c("cell","celltype","dataset","sample","stage")], by=c("cell"))
    
    # define color variable
    if (input$colourby == "gene_expression") {
      tmp <- data.table(
        cell = colnames(link_rna_expr),
        color = as.numeric(link_rna_expr[input$gene_umap,])
      )
      to.plot <- to.plot %>% merge(tmp,by="cell") %>% setorder(color)
    } else {
      to.plot$color <- to.plot[[input$colourby]]
    }

    ## Plot PAGA ##
    
    # celltypes <- sapply(paga$val,"[[","vertex.names")
    # alphas <- rep(1.0,length(celltypes)); names(alphas) <- celltypes
    # sizes <- rep(8,length(celltypes)); names(sizes) <- celltypes
    # 
    # p.paga <- ggnet2_interactive(
    #   net = paga,
    #   mode = c("x", "y"),
    #   color = celltype_colours[celltypes],
    #   node.alpha = alphas,
    #   node.size = sizes,    
    #   edge.size = 0.15,
    #   edge.color = "grey",
    #   label = TRUE,
    #   label.size = 3
    # )
    
    ## Plot UMAP ##
    
    # Subset UMAP for faster visualisation
    if (input$subset_cells_umap) {
      to.plot.umap <- to.plot[sample(1:.N, size=.N/3)]      
    } else {
      to.plot.umap <- to.plot
    }
    
    p.umap <- ggplot(to.plot.umap, aes(x = UMAP1, y = UMAP2, color = color)) +
      # geom_point(size = 1, alpha = 0.9) +
      geom_point_interactive(aes(tooltip = celltype, data_id = celltype), size=0.75, alpha=0.9) +
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
          axis.title = element_text(size = 10, color="black"),
          axis.text.y = element_text(size = 9, color="black"),
          axis.text.x = element_text(size = 9, color="black", angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }

    layout <- "
    BBB
    BBB
    BBB
    BBB
    CCC
    "
    girafe(
      # code = print(p.paga+p.umap+p2 + plot_layout(design = layout)),
      code = print(p.umap+p2 + plot_layout(design = layout)),
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
    # shiny::validate(need(input$gene_umap%in%genes, "" ))
    plot_UMAP()
  })
 

  ##################################
  ## Gene expression (pseudobulk) ##
  ##################################
  
  updateSelectizeInput(session = session, inputId = 'gene_pseudobulk', choices = genes, server = TRUE, selected = "T") 
  
  # TO-DO:
  # - fix order of classes 
  # - fix colors
  # - fix facet text
  # - add option to add number of cells for each barplot
  # - enable single-cell vs pseudobulk option
  
  plot_gene_expr_pseudobulk <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$gene_pseudobulk <- "T"
    # input$dataset_gene_expr_pseudoubulk <- c("KO", "CRISPR")
    # input$classes_gene_expr_pseudobulk <- classes
    # input$celltypes_gene_expr_pseudobulk <- celltypes
    ## END TEST ##
    
    sce.pseudobulk <- sce.pseudobulk[,sce.pseudobulk$class%in%input$classes_gene_expr_pseudobulk]
    sce.pseudobulk <- sce.pseudobulk[,sce.pseudobulk$celltype%in%input$celltypes_gene_expr_pseudobulk]
    sce.pseudobulk <- sce.pseudobulk[,sce.pseudobulk$dataset%in%input$dataset_gene_expr_pseudoubulk]
    
    ##############
    ## Barplots ##
    ##############
    
    to.plot <- data.table(
      # sample = colnames(sce.pseudobulk),
      expr = logcounts(sce.pseudobulk)[input$gene_pseudobulk,],
      class = sce.pseudobulk$class,
      sample = factor(sce.pseudobulk$sample, levels=names(class_colors)),
      celltype = factor(sce.pseudobulk$celltype, levels=input$celltypes_gene_expr_pseudobulk)
      # dataset = sce.pseudobulk$dataset
    )
    
    # %>% .[,celltype:=factor(celltype,levels=celltypes)] %>% 
      # .[,class:=factor(class,levels=opts$classes)]
    
    to.plot[expr==0,expr:=0.10]
    
    to.plot.means <- to.plot[,.(expr=mean(expr),sd=sd(expr)), by=c("class","celltype")]
    
    p.barplots <- ggplot(to.plot.means, aes(x=class, y=expr, fill=class)) +
      geom_bar(stat="identity", color="black", width=0.80) +
      geom_jitter(size=1, alpha=0.50, width=0.15, shape=21, data=to.plot) +
      geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.25, alpha=0.75, size=0.5) +
      facet_wrap(~celltype, scales="fixed") +
      scale_fill_manual(values=class_colors) +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",input$gene_pseudobulk)) +
      guides(x = guide_axis(angle = 90)) +
      theme(
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=rel(0.9)),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "none"
        # legend.title = element_blank(),
        # legend.text = element_text(size=rel(0.85))
      )
    
    if (input$add_number_observations_pseudobulk) {
      tmp <- max(to.plot$expr)+0.05
      p.barplots <- p.barplots +
        stat_summary(fun.data = function(x){ return(c(y = tmp, label = length(x))) }, geom = "text", size=2.75, data=to.plot)
    }
    
  # p.barplots <- cowplot::plot_grid(plotlist=p_list, ncol = 2)
  
  return(p.barplots)
})
  
    
  output$plot_gene_expr_pseudoubulk = renderPlot({
    shiny::validate(need(input$gene_pseudobulk%in%genes, "" ))
    plot_gene_expr_pseudobulk()
  })
    
    
  #############################
  ## Differential expression ##
  #############################
  
  
  # Update feature upon click
  observeEvent(input$diff_plot_selected, {
    updateSelectizeInput(session = session, inputId = 'diff_feature', choices = genes, server = TRUE, selected = input$diff_plot_selected)
  })
  
  plot_differential <- reactive({
    
    ## START TEST ##
    input <- list()
    input$diff_class <- "Dnmt1_KO"
    input$diff_celltype <- "Gut"
    input$diff_gene <- "Foxa2"
    input$diff_resolution <- "Cells"
    input$diff_range <- c(-8,8)
    input$diff_min_log_pval <- 25
    ## END TEST ##
    
    ## Volcano ##
    
    # Fetch data
    if (input$diff_resolution=="Cells") {
      diff.dt <- fread(file.path(data_folder,sprintf("differential/cells/%s/%s_WT_vs_%s.txt.gz",input$diff_class,input$diff_celltype,input$diff_class))) %>%
        .[,c("groupA_N","groupB_N"):=NULL]
    } else if (input$diff_resolution=="Pseudobulk") {
      stop("Not implemented")
      # diff.dt <- fread(file.path(data_folder,sprintf("differential/pseudobulk/%s/%s_WT_vs_%s.txt.gz",input$diff_class,input$diff_celltype,input$diff_class)))
    }
    # diff.dt <- diff.dt[gene%in%genes]
    
    if (!input$diff_gene%in%diff.dt$gene) {
      diff.dt <- rbind(diff.dt, data.table(gene=input$diff_gene, logFC=0, padj_fdr=1))
    }
    
    to.plot <- diff.dt %>% 
      .[logFC>=input$diff_range[1] & logFC<=input$diff_range[2]] %>% 
      .[is.na(logFC),logFC:=0] %>%
      .[is.na(padj_fdr),padj_fdr:=1] %>%
      .[,log_pval:=-log10(padj_fdr+1e-150)]
    
    
    # xlim_min <- min(to.plot$logFC); xlim_max <- max(to.plot$logFC)
    xlim_min <- -8; xlim_max <- 8; margin_text <- 2; margin_dots <- 0.15
    to.plot[logFC<=xlim_min,logFC:=xlim_min]; to.plot[logFC>=xlim_max,logFC:=xlim_max]
    
    # negative_hits <- to.plot[sig==TRUE & r<0 & r>input$rna_vs_chromvar_cor_range[1],gene]
    # positive_hits <- to.plot[sig==TRUE & r>0 & r<input$rna_vs_chromvar_cor_range[2],gene]
    
    # to.plot.subset <- rbind(
    #   to.plot[log_pval>=10],
    #   to.plot[log_pval<=10][sample.int(.N,5e3)]
    # )
    to.plot.subset <- to.plot[log_pval>=input$diff_min_log_pval]
    if (!input$diff_gene%in%to.plot.subset$gene) {
      to.plot.subset <- rbind(to.plot.subset,to.plot[gene==input$diff_gene])
    }
    
    ylim_min <- min(to.plot.subset$log_pval,na.rm=T); ylim_max <- max(to.plot$log_pval,na.rm=T)
    
    p.volcano <- ggplot(to.plot.subset, aes(x=logFC, y=log_pval)) +
      geom_segment(x=0, xend=0, y=input$diff_min_log_pval, yend=ylim_max, color="orange", size=0.25) +
      geom_jitter_interactive(aes(fill=log_pval, tooltip=gene, data_id=gene, alpha=log_pval, size=log_pval, onclick=gene), width=0.03, height=0.03, shape=21) + 
      geom_point_interactive(aes(tooltip=gene, data_id=gene, onclick=gene), size=7, fill="green", shape=21, data=to.plot.subset[gene==input$diff_gene]) + 
      geom_hline(yintercept=input$diff_min_log_pval, linetype="dashed") +
      ggrepel::geom_text_repel(aes(label=gene), size=5, data=to.plot[gene==input$diff_gene]) +
      scale_fill_gradient(low = "gray80", high = "red") +
      scale_alpha_continuous(range=c(0.25,1)) +
      scale_size_continuous(range=c(0.15,3.5)) +
      scale_x_continuous(limits=c(xlim_min-margin_dots,xlim_max+margin_dots)) +
      scale_y_continuous(limits=c(ylim_min,ylim_max+3)) +
      annotate("text", x=0, y=ylim_max+3, size=5, label=sprintf("(%d)", nrow(to.plot.subset))) +
      # coord_cartesian(xlim=c(xlim_min,xlim_max)) +
      # annotate("text", x=xlim_min-0.05, y=ylim_max+2, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      # annotate("text", x=xlim_max+0.05, y=ylim_max+2, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      annotate("text", x=xlim_min+margin_text, y=1, size=5, label=sprintf("Higher in %s",input$diff_celltypeA)) +
      annotate("text", x=xlim_max-margin_text, y=1, size=5, label=sprintf("Higher in %s",input$diff_celltypeB)) +
      labs(x=ifelse(input$diff_modality=="RNA","Differential expression","Differential accessibility"), y=expression(paste("-log"[10],"(p.value)"))) +
      theme_classic() +
      theme(
        axis.text = element_text(size=rel(1.25), color='black'),
        axis.title = element_text(size=rel(1.25), color='black'),
        legend.position="none"
      )
    
    ## Plot expression/accessibility values ##
    
    # Fetch data
    stop()
    if (input$diff_resolution=="Cells") {
      expr.dt <- data.table(
        # cell = 
        expr = as.numeric(rna_expr_cells.array[input$diff_gene,])
        # celltype = cell_metadata.dt[colnames(rna_expr_cells.array),celltype]
      )
    } else if (input$diff_resolution=="Pseudobulk") {
      expr.dt <- data.table(
        expr = as.numeric(rna_expr_pseudobulk_replicates.mtx[input$diff_gene,]),
        sample = colnames(rna_expr_pseudobulk_replicates.mtx)
      ) %>% .[,celltype:=strsplit(sample,"-") %>% map_chr(1)]
    }
    
    expr.dt <- expr.dt %>% 
      .[celltype%in%c(input$diff_celltypeA,input$diff_celltypeB)] %>% 
      .[,celltype:=factor(celltype,levels=c(input$diff_celltypeA,input$diff_celltypeB))]
    
    
    if (input$diff_resolution=="Cells") {
      
      p.expr <- ggplot(expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
        geom_jitter(size=2, width=0.05, alpha=0.5, shape=21) +
        geom_violin(scale="width", alpha=0.40) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.70) +
        stat_summary(fun.data = function(x) { return(c(y = max(expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=5) +
        scale_fill_manual(values=celltype_colours[c(input$diff_celltypeA,input$diff_celltypeB)]) +
        labs(x="", y=sprintf("%s expression",input$diff_gene), title=input$diff_gene) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.25)),
          axis.text.x = element_text(colour="black",size=rel(1.25)),
          axis.text.y = element_text(colour="black",size=rel(1.25)),
          axis.title.y = element_text(colour="black",size=rel(1.25)),
          axis.ticks.x = element_blank(),
          legend.position = "none"
        )
      
    } else if (input$diff_resolution=="Pseudobulk") {
      
      tmp <- expr.dt[,.(expr=mean(expr), sd=sd(expr)), by="celltype"]
      
      p.expr <- ggplot(expr.dt, aes(x=celltype, y=expr, fill=celltype)) +
        geom_bar(stat="identity", color="black", alpha=0.9, data=tmp) +
        geom_jitter(size=2.5, alpha=0.9, width=0.08, shape=21) +
        geom_errorbar(aes(ymin=expr-sd, ymax=expr+sd), width=0.25, alpha=1, size=0.6, data=tmp) +
        stat_summary(fun.data = function(x) { return(c(y = max(expr.dt$expr)+0.5, label = length(x)))}, geom = "text", size=5) +
        scale_fill_manual(values=celltype_colours[c(input$diff_celltypeA,input$diff_celltypeB)]) +
        labs(x="",y=sprintf("%s expression",input$diff_gene), title=input$diff_gene) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust=0.5, size=rel(1.25)),
          axis.text.x = element_text(colour="black",size=rel(1.25)),
          axis.text.y = element_text(colour="black",size=rel(1.25)),
          axis.title.y = element_text(colour="black",size=rel(1.25)),
          axis.ticks.x = element_blank(),
          legend.position = "none"
        )
    }
    
    # girafe call
    girafe(
      code = print(p.volcano + p.expr + plot_layout(ncol=2, widths=c(2,1))),
      width_svg = 14, height_svg = 9,
      options = list(
        opts_zoom(min = 1, max = 5),
        opts_sizing(rescale = TRUE),
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer")
      )
    ) %>% return(.)
    
  })
  
  output$diff_plot <- renderGirafe({
    shiny::validate(need(input$diff_celltypeA%in%celltypes,"Please select celltype A"))
    shiny::validate(need(input$diff_celltypeB%in%celltypes,"Please select celltype B"))
    shiny::validate(need(input$diff_celltypeA!=input$diff_celltypeB,"Celltype A and B must be different"))
    shiny::validate(need(input$diff_resolution%in%c("Cells","Pseudobulk"),"Please select data resolution (cell, pseudobulk)"))
    shiny::validate(need(input$diff_feature%in%genes,"Please select a gene from our annotation"))
    plot_differential()
  })
  

  
  plot_diff_heatmap <- reactive({
    
    # TO-DO:
    # - split by dataset
    # - click on a gene allows you to see the barplots/boxplots
    
    
    ## START TEST ##
    # input <- list()
    # input$genes_diff_heatmap <- c("T","Sox2","Foxa2","Hoxd9")
    # input$celltypes_diff_heatmap <- celltypes_pseudobulk
    # input$classes_diff_heatmap <- classes[classes!="WT"]
    ## END TEST ##
    
    to.plot <- expand.grid(input$celltypes_diff_heatmap, input$classes_diff_heatmap, input$genes_diff_heatmap) %>% 
      as.data.table %>% setnames(c("celltype","class","gene")) %>%
      merge(diff_pseudobulk.dt[gene%in%input$genes_diff_heatmap & celltype%in%input$celltypes_diff_heatmap & class%in%input$classes_diff_heatmap] , by=c("celltype","class","gene"), all.x=T)
    
    # define order of plotting
    to.plot %>% .[,class:=factor(class,levels=rev(classes[classes%in%input$classes_diff_heatmap]))]
    to.plot %>% .[,celltype:=factor(celltype,levels=celltypes[celltypes%in%input$celltypes_diff_heatmap])]
    
    # filter entries with lots of NAs
    to.plot <- to.plot[,foo:=mean(is.na(diff)),by=c("celltype","gene")] %>% .[foo<1] %>% .[,foo:=NULL]
    
    p <- ggplot(to.plot, aes(x=celltype, y=class, fill=diff)) +
      # geom_tile(color="black") +
      geom_tile_interactive(aes(tooltip=diff), color="black") +
      facet_wrap(~gene) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = 'gray70' ) +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        axis.text = element_text(color="black", size=rel(0.75)),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black", size=rel(1.25)),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank()
      )
    
    girafe(
      code = print(p),
      width_svg = 13, height_svg = 9,
      options = list( 
        opts_sizing(rescale = FALSE)
      )
    ) %>% return(.)
    
  })
  
  output$plot_diff_heatmap = renderGirafe({
    plot_diff_heatmap()
  })
  
  
  ##########################
  ## Celltype proportions ##
  ##########################
  
  plot_celltype_proportions <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$class_celltype_proportions <- "WT"
    # input$dataset_celltype_proportions <- c("KO", "CRISPR")
    # input$visualisation_type_celltype_proportions <- "barplots"
    # input$remove_extraembryonic <- FALSE
    ## END TEST ##
    
    to.plot <- cell_metadata.dt %>% 
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
    
    # Plot
    if (input$visualisation_type_celltype_proportions=="Barplots per sample") {
      
      # Define sample order
      to.plot[,sample:=sprintf("(%s) %s",dataset,sample)]
      
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
          legend.position = "none",
          # legend.key.size = unit(0.50, "cm"),
          # strip.background = element_blank(),
          # strip.text = element_text(color="black", size=rel(0.9)),
          axis.title.x = element_text(color="black", size=rel(0.75)),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=rel(0.5), color="black"),
          axis.text.x = element_text(size=rel(1), color="black")
        )
    } else if (input$visualisation_type_celltype_proportions=="Boxplots per class") {
      
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
    # shiny::validate(need(input$gene_umap%in%genes, "" ))
    plot_celltype_proportions()
  })
  
  
  ######################################
  ## Celltype proportions comparisons ##
  ######################################
  
  plot_celltype_proportions_comparisons <- reactive({
    
    ## START TEST ##
    # input <- list()
    # input$class_celltype_comparisons <- c("Dnmt1_KO","Dnmt3a_KO")
    # input$dataset_celltype_comparisons <- c("KO", "CRISPR")
    # input$remove_extraembryonic_celltype_comparisons <- TRUE
    ## END TEST ##
    
    # remove small embryos
    # if (opts$remove.small.embryos) {
    #   opts$min.cells <- 1500
    #   cell_metadata.dt <- cell_metadata.dt %>%
    #     .[,N:=.N,by="sample"] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
    # }
    
    # Select KO samples    
    cell_metadata_filt.dt <- cell_metadata.dt[dataset%in%input$dataset_celltype_comparisons]
    
    if (input$remove_extraembryonic_celltype_comparisons) {
      cell_metadata_filt.dt <- cell_metadata_filt.dt[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
    }
    
    # calculate celltype proportions for WT samples
    wt_proportions.dt <- cell_metadata_filt.dt %>%
      .[class=="WT"] %>%
      .[,ncells:=.N] %>%
      .[,.(proportion=.N/unique(ncells), N=round(.N/length(unique(sample)))), by="celltype"]
    
    # calculate celltype proportions for KO samples
    ko_proportions_per_sample.dt <- cell_metadata_filt.dt %>%
      .[class%in%input$class_celltype_comparisons] %>%
      setkey(celltype,sample) %>%
      .[CJ(celltype,sample, unique = TRUE), .N, by = .EACHI] %>%
      merge(unique(cell_metadata_filt.dt[,c("sample","class")]), by="sample") %>%
      .[,ncells:=sum(N), by="sample"] %>% .[,proportion:=(N+1)/ncells]
    
    # Calculate difference in celltype proportions between KO and WT
    proportions_per_sample.dt <- merge(
      ko_proportions_per_sample.dt, 
      wt_proportions.dt, 
      by = c("celltype"), allow.cartesian=T, suffixes = c(".ko",".wt")
    ) %>% .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] 
    
    ylimits <- max(abs(proportions_per_sample.dt[!is.infinite(diff_proportion),diff_proportion]))
    
    # Filter cell types with at least N cells
    celltypes.to.plot <- proportions_per_sample.dt %>%
      .[,.(N=sum(N.ko)+sum(N.wt)),by=c("class","celltype")] %>% 
      .[N>=50,celltype] %>% as.character
    
    to.plot <- proportions_per_sample.dt %>%
      .[celltype%in%celltypes.to.plot] %>%
      merge(unique(cell_metadata_filt.dt[,c("sample","dataset")]),by="sample")
    
    # Define cell type order
    # if (length(input$class_celltype_comparisons))>=
    celltype.order <- to.plot %>%
      .[,mean(diff_proportion),by="celltype"] %>% setorder(-V1) %>% .$celltype
    to.plot <- to.plot %>% .[,celltype:=factor(celltype,levels=celltype.order)]
    
    # Plot
    p <- ggplot(to.plot, aes(x=celltype, y=diff_proportion)) +
      geom_point(aes(fill = celltype), shape=21, size=1.5) +
      geom_boxplot_interactive(aes(tooltip=celltype, fill=celltype), alpha=0.5) +
      coord_flip(ylim=c(-ylimits,ylimits)) +
      geom_hline(yintercept=0, linetype="dashed", size=0.5) +
      scale_fill_manual(values=celltype_colours, drop=F) +
      theme_classic() +
      labs(y="Difference in proportions (log2)", x="") +
      theme(
        legend.position = "none",
        strip.text = element_text(color="black", size=rel(1.25)),
        # axis.title = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black")
      )
    
    # Split by data set
    if (input$split_by_dataset_celltype_comparisons) {
      p <- p + facet_wrap(~dataset,nrow=1)
    }
      
    girafe(
      code = print(p),
      width_svg = 11, height_svg = 9,
      options = list(
        opts_sizing(rescale = FALSE)
      )
    ) %>% return(.)
    
  })
  
  output$plot_celltype_proportions_comparisons = renderGirafe({
    plot_celltype_proportions_comparisons()
  })

}
