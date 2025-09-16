#----------------------------------#
#-- Fit C-side to simulated data --#
#----------------------------------#

data_dir <- "~MERFISH/"
library(tidyverse)
library(Seurat)
library(peakRAM)
library(spacexr)
library(Matrix)
library(doParallel)
benchmark_and_return <- function(expr) {
  expr <- substitute(expr) 
  calling_env <- parent.frame()       
  eval_env <- new.env(parent = calling_env) 
  mem <- peakRAM({
    eval_env$.__res__ <- eval(expr, envir = eval_env)
  })
  list(result = eval_env$.__res__, memory = mem)
}

load( paste0(data_dir, "Allen/WholeBrain_seurat_distance.Rdata"))

# Getting reference with other cell types
load(paste0(data_dir, "Allen/WholeBrain_seurat.Rdata"))
load(file.path(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
colnames(seurat_allen@meta.data)[4] <- "cell_ID"

# Filtering gene
filteredgenes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1) %>%pull(target)
seurat_allen <- subset(seurat_allen, features = filteredgenes)
seurat_allen <- subset(seurat_allen, subset = division == "Isocortex")
seurat_allen <- subset(seurat_allen, subset = class %in% c("01 IT-ET Glut","NP-CT-L6b Glut","30 Astro-Epen ","33 Vascular","34 Immune"))
# add DE variable
seurat_allen@meta.data <- seurat_allen@meta.data %>% 
  left_join(dataDE@meta.data[,c("cell_ID", "kernelMicroglia1000")])

# Deconvolution step
celltypes <- factor(seurat_allen@meta.data$class)
nUMI <- seurat_allen@meta.data$nCount_RNA
names(celltypes) <- names(nUMI) <- seurat_allen@meta.data$cell_ID
reference <- spacexr::Reference(seurat_allen@assays$RNA$counts, celltypes, nUMI = nUMI, min_UMI = 0)

coords <- seurat_allen@meta.data[,c("x","y")]
rownames(coords) <- seurat_allen@meta.data$cell_ID
spaceRNA <- SpatialRNA(coords = coords
                       ,counts = seurat_allen@assays$RNA$counts
                       ,nUMI = nUMI)

myRCTD <- create.RCTD(spaceRNA, reference,max_cores = 5)
myRCTD <-benchmark_and_return(run.RCTD(myRCTD, doublet_mode = 'doublet'))

X <- seurat_allen@meta.data$kernelMicroglia1000
names(X) <-  seurat_allen@meta.data$cell_ID
X[is.na(X)]<-0
results_cside <- benchmark_and_return(run.CSIDE.single(
  myRCTD$result,
  explanatory.variable = X ,
  cell_types = "01 IT-ET Glut",
  cell_type_threshold = 0,
  gene_threshold = 0,
  doublet_mode = T,
  log_fc_thresh = 0,
  weight_threshold = 0,
  fdr=1,
  medv = 0
))
save(myRCTD,results_cside,
     file = paste0(data_dir, "ResultsDistance/ResultsCside.Rdata"))
