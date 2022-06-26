
#####################
## Define settings ##
#####################

if (Sys.info()[['nodename']]=="BI2404M") {
  data_folder <- "/Users/argelagr/shiny_dnmt/data"
} else if (Sys.info()[['nodename']]=="Ricards-MacBook-Pro.local") {
  data_folder <- "/Users/rargelaguet/shiny_dnmt/data"
}

###########################
## Load global variables ##
###########################

# Load genes
genes <- fread(paste0(data_folder,"/genes.txt"), header=F)[[1]]
cells <- fread(paste0(data_folder,"/cells_rna.txt"), header=F)[[1]]

########################
## Load cell metadata ##
########################

sample_metadata <- fread(paste0(data_folder,"/cell_metadata.txt.gz")) %>% 
  .[,celltype := factor(celltype, levels = celltypes, ordered = TRUE)]
# .[,sample = factor(sample, levels = names(samples), ordered = TRUE)

#########################
## Load RNA expression ##
#########################

# Load gene expression matrix
link_rna_expr = HDF5Array(file = paste0(data_folder,"/rna_expr.hdf5"), name = "expr_logcounts")
colnames(link_rna_expr) <- cells
rownames(link_rna_expr) <- genes

# Load pseudobulk SingleCellExperiment object
# sce.pseudobulk <- readRDS(file.path(data_folder,"SingleCellExperiment_pseudobulk_class_celltype_dataset.rds"))
sce.pseudobulk <- readRDS(file.path(data_folder,"SingleCellExperiment_shiny.rds"))
celltypes_pseudobulk <- celltypes[celltypes%in%unique(sce.pseudobulk$celltype)]

#####################
## Load PAGA graph ##
#####################

# paga <- readRDS(paste0(data_folder,"/paga_network.rds"))

connectivity.mtx <- fread(file.path(data_folder,"paga/paga_connectivity.csv")) %>%
  matrix.please %>% .[celltypes,celltypes]

coordinates.mtx <- fread(file.path(data_folder,"paga/paga_coordinates.csv")) %>% 
  matrix.please %>% .[celltypes,]

# Parse data
connectivity.mtx[connectivity.mtx<0.20] <- 0
connectivity.mtx[connectivity.mtx>=0.20] <- 1

# Create igraph object
igraph.paga <- graph_from_adjacency_matrix(connectivity.mtx, mode = "undirected")

# Create tbl_graph object
igraph.paga.tbl <- as_tbl_graph(igraph.paga) %>%
  activate(nodes) %>%
  mutate(celltype=rownames(connectivity.mtx)) %>%
  mutate(x=coordinates.mtx[,1]) %>% mutate(y=coordinates.mtx[,2])

# Create network object
net.paga = network(connectivity.mtx)
net.paga %v% "x" = coordinates.mtx[,1]
net.paga %v% "y" = coordinates.mtx[,2]

######################################
## Load differential RNA expression ##
######################################

diff_pseudobulk.dt <- fread(file.path(data_folder,"diff_pseudobulk.txt.gz")) %>% .[celltype%in%celltypes_pseudobulk]

###################################
## Load dimensionality reduction ##
###################################

# Load UMAPs
umap_list <- classes %>% map(function(i) {
  fread(file.path(data_folder,sprintf("/dimensionality_reduction/%s/umap.txt.gz",i)))
}) %>% set_names(classes)

#########################
## Load reference UMAP ##
#########################

# umap_reference.dt <- fread(paste0(data_folder,"/mapping/umap_coordinates.txt.gz"))