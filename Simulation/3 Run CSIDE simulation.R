#----------------------------------#
#-- Fit C-side to simulated data --#
#----------------------------------#

dir.out = "~/Simulation/"
library(tidyverse)
library(Seurat)
library(peakRAM)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
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

load(paste0(dir.out,paste0("SimulatedY.Rdata")))
n <- nrow(data)
data$cell_ID <- as.character(data$cell_ID)
data$niche <- data$nichenum
data$nichenum <- ifelse(data$niche == "stroma",1,0)
rownames(Ycounts) <- data$cell_ID

# Generating another cell type
load(paste0(dir.out,"Lung5-3data.Rdata"))
data2 <-  fulldata %>% mutate(auxID = 1:nrow(fulldata)) %>%
  filter(cell_type == "plasmablast",
         niche %in% c("myeloid-enriched stroma","stroma"))
data2$nichenum <- ifelse(data2$niche == "stroma",1,0)
n2 <- nrow(data2)
p<-nrow(params)
Ycounts2 <- do.call(cbind, lapply(1:p, function(gene){
  set.seed(31*gene)
  rpois(n2, 10)
}))
colnames(Ycounts2) <- params$gene_name
rownames(Ycounts2) <- as.character(data2$cell_ID)

# Combining with macrophages
data$clust <- "main"
data2$clust <- "aux"
Ycounts <- rbind(Ycounts, Ycounts2)

data <- rbind(data, data2)
rownames(data) <- data$cell_ID

# Get reference dataset
celltypes <- factor(data$clust)
nUMI <- rowSums(Ycounts)
names(celltypes) <- names(nUMI) <- data$cell_ID
reference <- spacexr::Reference(t(Ycounts), celltypes, nUMI = nUMI, min_UMI = 0)

coords <- data[,c("sdimx", "sdimy")]
rownames(coords) <- data$cell_ID
spaceRNA <- SpatialRNA(coords = coords,
                       counts = t(Ycounts),
                       nUMI = nUMI)

myRCTD <- create.RCTD(spaceRNA, reference, gene_cutoff = 0,
                      fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg=0,
                      UMI_min = 0, counts_MIN = 0, CELL_MIN_INSTANCE =0)
myRCTD <-benchmark_and_return(run.RCTD(myRCTD, doublet_mode = 'doublet'))

# Run C-side
X <- data$nichenum
names(X) <- data$cell_ID
results_cside <- benchmark_and_return(run.CSIDE.single(
  myRCTD$result,
  explanatory.variable = X ,
  cell_types = "main",
  cell_type_threshold = 0,
  gene_threshold = 0,
  doublet_mode = T,
  log_fc_thresh = 0,
  weight_threshold = 0,
  fdr=1,
  medv = 0
))
save(myRCTD,results_cside,
     file = paste0(dir.out, "Results/PartialSimulation/ResultsCside.Rdata"))

