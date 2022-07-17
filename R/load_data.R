
#####################
## Define settings ##
#####################

# if (Sys.info()[['nodename']]=="BI2404M") {
#   data_folder <- "/Users/argelagr/shiny_dnmt/data"
# } else if (Sys.info()[['nodename']]=="rargelaguet.local") {
#   data_folder <- "/Users/rargelaguet/shiny_dnmt/data"
# }

###########################
## Load global variables ##
###########################

genes <- fread(paste0(data_folder,"/rna_expression/genes.txt"), header=F)[[1]]
cells <- fread(paste0(data_folder,"/rna_expression/rna_expr_cells.txt"), header=F)[[1]]

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(paste0(data_folder,"/cell_metadata.txt.gz")) %>% 
  # .[celltype%in%celltypes] %>%
  .[,celltype := factor(celltype, levels = celltypes, ordered = TRUE)] %>%
  setkey(cell)
# .[,sample = factor(sample, levels = names(samples), ordered = TRUE)

# unique(cell_metadata.dt$celltype)[!unique(cell_metadata.dt$celltype)%in%celltypes]
# sum(is.na(cell_metadata.dt$celltype))

#################################
## Load RNA expression (cells) ##
#################################

rna_expr_cells.array = HDF5Array(file.path(data_folder,"rna_expression/rna_expr_cells.hdf5"), name = "rna_expr_logcounts")
colnames(rna_expr_cells.array) <- cells
rownames(rna_expr_cells.array) <- genes

stopifnot(cell_metadata.dt$cell==colnames(rna_expr_cells.array))

######################################
## Load RNA expression (pseudobulk) ##
######################################

sample_metadata_pseudobulk_list <- list()
rna_expr_pseudobulk_list <- list()

# load "class_sample_celltype_dataset"
sample_metadata_pseudobulk_list[["class_sample_celltype_dataset"]] <- fread(file.path(data_folder,"rna_expression/pseudobulk/class_sample_celltype_dataset/sample_metadata.txt.gz"))
rna_expr_pseudobulk_list[["class_sample_celltype_dataset"]] <- readRDS(file.path(data_folder,"rna_expression/pseudobulk/class_sample_celltype_dataset/rna_expr.rds"))

# load "class_sample_celltype"
sample_metadata_pseudobulk_list[["class_sample_celltype"]] <- fread(file.path(data_folder,"rna_expression/pseudobulk/class_sample_celltype/sample_metadata.txt.gz"))
rna_expr_pseudobulk_list[["class_sample_celltype"]] <- readRDS(file.path(data_folder,"rna_expression/pseudobulk/class_sample_celltype/rna_expr.rds"))

# sce.pseudobulk <- readRDS(file.path(data_folder,"rna_expression/SingleCellExperiment_pseudobulk.rds"))
# sce.pseudobulk <- readRDS(file.path(data_folder,"rna_expression/SingleCellExperiment_pseudobulk.rds"))
# celltypes_pseudobulk <- celltypes[celltypes%in%unique(sce.pseudobulk$celltype)]

#####################
## Load PAGA graph ##
#####################

# # paga <- readRDS(paste0(data_folder,"/paga_network.rds"))
# 
# connectivity.mtx <- fread(file.path(data_folder,"paga/paga_connectivity.csv")) %>%
#   matrix.please %>% .[celltypes,celltypes]
# 
# coordinates.mtx <- fread(file.path(data_folder,"paga/paga_coordinates.csv")) %>% 
#   matrix.please %>% .[celltypes,]
# 
# # Parse data
# connectivity.mtx[connectivity.mtx<0.20] <- 0
# connectivity.mtx[connectivity.mtx>=0.20] <- 1
# 
# # Create igraph object
# igraph.paga <- graph_from_adjacency_matrix(connectivity.mtx, mode = "undirected")
# 
# # Create tbl_graph object
# igraph.paga.tbl <- as_tbl_graph(igraph.paga) %>%
#   activate(nodes) %>%
#   mutate(celltype=rownames(connectivity.mtx)) %>%
#   mutate(x=coordinates.mtx[,1]) %>% mutate(y=coordinates.mtx[,2])
# 
# # Create network object
# net.paga = network(connectivity.mtx)
# net.paga %v% "x" = coordinates.mtx[,1]
# net.paga %v% "y" = coordinates.mtx[,2]


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

umap_reference.dt <- fread(paste0(data_folder,"/mapping/umap_coordinates.txt.gz"))

#########################
## Repetitive elements ##
#########################

repeats_expr.dt <- fread(file.path(data_folder,"repeats/repeats_expr.txt.gz"))
repeats_diff_expr.dt <- fread(file.path(data_folder,"repeats/repeats_diff_expr.txt.gz"))
