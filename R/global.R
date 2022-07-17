############################################
## Define computer (for testing purposes) ##
############################################

if (Sys.info()[['nodename']]=="BI2404M") {
  data_folder <- "/Users/argelagr/shiny_dnmt/data"
} else if (grepl("rargelaguet",Sys.info()[['nodename']])) {
  # data_folder <- "/Users/rargelaguet/shiny_dnmt/data"
  data_folder <- "/Users/rargelaguet/data/10x_gastrulation_DNMTs/shiny"
}

###############
## libraries ##
###############

library(R.utils)
library(HDF5Array)
library(data.table)
library(purrr)
library(DT)

# shiny
library(shiny)
library(shinyFiles)
library(shinythemes)
library(ggiraph)

# general viz
library(GGally)
library(cowplot)
library(ggrepel)
library(ggplot2)
require(patchwork)
require(ggpubr) # to do: remove this dependency?

# graph viz
# require(visNetwork)
# library(sna)
library(network)
library(ggraph)
library(igraph)
library(tidygraph)

######################
## Global variables ##
######################

classes <- c(
  "WT", 
  "Dnmt3a_KO", 
  "Dnmt3b_KO",
  "Dnmt1_KO"
)

celltypes <- c(
  # "Epiblast",
  # "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  # "Nascent_mesoderm",
  # "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  # "Erythroid3",
  "Erythroid",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm"
  # "Parietal_endoderm"
)

samples <- c("Dnmt3a_KO_1", "Dnmt3a_KO_2", "WT_1", "WT_2", "WT_3", "WT_4", 
             "WT_5", "WT_6", "WT_7", "Dnmt3a_KO_13", "Dnmt3a_KO_14", "Dnmt3b_KO_11", 
             "Dnmt3b_KO_1", "Dnmt3b_KO_2", "Dnmt3ab_KO_1", "Dnmt3b_KO_12", 
             "Dnmt1_KO_1", "Dnmt1_KO_2", "Dnmt1_KO_3", "Dnmt1_KO_4", "Dnmt1_KO_5", 
             "Dnmt1_KO_6", "Dnmt1_KO_7", "Dnmt1_KO_8", "Dnmt1_KO_9", "Dnmt1_KO_10", 
             "Dnmt1_KO_11", "Dnmt1_KO_12", "Dnmt1_KO_13", "Dnmt1_KO_14", "Dnmt1_KO_15", 
             "Dnmt3a_KO_3", "Dnmt3a_KO_4", "Dnmt3a_KO_5", "Dnmt3a_KO_6", "Dnmt3a_KO_7", 
             "Dnmt3a_KO_8", "Dnmt3a_KO_9", "Dnmt3a_KO_10", "Dnmt3a_KO_11", 
             "Dnmt3a_KO_12", "Dnmt3b_KO_3", "Dnmt3b_KO_4", "Dnmt3b_KO_5", 
             "Dnmt3b_KO_6", "Dnmt3b_KO_7", "Dnmt3b_KO_8", "Dnmt3b_KO_9", "Dnmt3b_KO_10", 
             "WT_8", "WT_9", "WT_10", "WT_11", "WT_12", "WT_13", "WT_14", 
             "WT_15", "WT_16", "WT_17")

repeat_classes <- c(
  "LINE_L1",
  "LINE_L2",
  "LTR_ERV1",
  "LTR_ERVK",
  "LTR_ERVL",
  "LTR_MaLR",
  "major_satellite",
  "minor_satellite",
  "rRNA",           
  "SINE_Alu_B1",
  "SINE_B2",
  "SINE_B4",
  "IAP"  
)

#####################
## Colour palettes ##
#####################

celltype_colours <- c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  # "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood_progenitors" = "#c9a997",
  # "Blood_progenitors_1" = "#f9decf",
  # "Blood_progenitors_2" = "#c9a997",
  "Erythroid" = "#EF4E22",
  # "Erythroid1" = "#C72228",
  # "Erythroid2" = "#f79083",
  # "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  # "Neurectoderm" = "#65A83E",
  "Rostral_neurectoderm" = "#65A83E",
  "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A"
)

celltype_palette = scale_color_manual(values = celltype_colours, name = "", drop=TRUE)
celltype_palette_fill = scale_fill_manual(values = celltype_colours, name = "", drop=TRUE)

sample_palette <- scale_color_brewer(palette="Paired")
sample_palette_fill <- scale_fill_brewer(palette="Paired")

dataset_palette <- scale_color_brewer(palette="Dark2")
dataset_palette_fill <- scale_fill_brewer(palette="Dark2")

rna_palette <- scale_color_gradient(low = "gray80", high = "red")

# stage_colours = c("E7.5" = "#FFFFBF", "E8.5" = "#3288BD")
# stage_palette = scale_color_manual(values = stage_colours, name = "stage")
# stage_palette_fill = scale_fill_manual(values = stage_colours, name = "stage")

class_colors <- c(
  "WT" = "#ffffb3", 
  "Dnmt3a_KO" = "#8dd3c7", 
  "Dnmt3b_KO" = "#fb8072", 
  "Dnmt1_KO" = "#80b1d3"
  # "Dnmt3ab_KO" = "#bebada"
)

###############
## Functions ##
###############

minmax.normalisation <- function(x) {
  return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}


matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
