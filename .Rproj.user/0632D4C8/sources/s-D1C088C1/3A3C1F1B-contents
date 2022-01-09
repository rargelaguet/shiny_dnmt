
# Load UMAPs
umap_list <- classes %>% map(function(i) {
  fread(file.path(basedir,sprintf("/dimensionality_reduction/%s/umap.txt.gz",i)))
}) %>% set_names(classes)


# Load PAGA
paga <- readRDS(paste0(basedir,"/paga_network.rds"))

# Load cell metadata
sample_metadata <- fread(paste0(basedir,"/cell_metadata.txt.gz")) %>% 
  setnames("celltype.mapped","celltype") %>% 
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[,celltype := factor(celltype, levels = celltypes, ordered = TRUE)]
  # .[,sample = factor(sample, levels = names(samples), ordered = TRUE)

# Load genes
genes <- fread(paste0(basedir,"/genes.txt"), header=F)[[1]]

# Load gene expression matrix
link_rna_expr = HDF5Array(file = paste0(basedir,"/rna_expr.hdf5"), name = "expr_logcounts")
colnames(link_rna_expr) <- fread(paste0(basedir,"/cells_rna.txt"), header=F)[[1]]
rownames(link_rna_expr) <- genes

