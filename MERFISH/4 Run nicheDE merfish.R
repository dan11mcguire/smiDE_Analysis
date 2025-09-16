#---------------------#
#- Run nicheDE Models -#
#---------------------#
#devtools::install_github("kaishumason/NicheDE") 
library(nicheDE)
library(tidyverse)
library(peakRAM)

benchmark_and_return <- function(expr) {
  expr <- substitute(expr)             
  calling_env <- parent.frame()        
  eval_env <- new.env(parent = calling_env) 
  mem <- peakRAM({
    eval_env$.__res__ <- eval(expr, envir = eval_env)
  })
  list(result = eval_env$.__res__, memory = mem)
}

data_dir <- "~MERFISH/"
ncores<- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')) -1
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
sigma <- 1000
dist <- c("negbin","gaus")[as.numeric(slurm_arrayid)]

# Load data
load(paste0(data_dir, "Allen/WholeBrain_seurat.Rdata"))
load(paste0(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
colnames(seurat_allen@meta.data)[4] <- "cell_ID"

# Filtering gene
filteredgenes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1) %>% pull(target)
seurat_allen <- subset(seurat_allen, features = filteredgenes)
seurat_allen <- subset(seurat_allen, subset = division == "Isocortex")

dataDE <- subset(seurat_allen, subset = class == "01 IT-ET Glut" | subclass == "334 Microglia NN")
dataDE@meta.data$class <- ifelse(dataDE@meta.data$class  == "01 IT-ET Glut", "Neurons","Microlglia")

# Create deconvolution matrix (dummy)
dummy_class <- model.matrix(~class-1, data = dataDE@meta.data)
colnames(dummy_class) <- gsub("class","",colnames(dummy_class))


if(dist == "negbin"){
  # Get average expression profile
  avg_expr <- CreateLibraryMatrix(t(dataDE@assays$RNA$counts),
                                  dataDE@meta.data[,c("cell_ID","class")])
  # NicheDE object
  NDE_obj <- CreateNicheDEObject(counts_mat = t(dataDE@assays$RNA$counts),
                                 coordinate_mat = dataDE@meta.data[,c("x","y")],
                                 library_mat = avg_expr[colnames(dummy_class),],
                                 deconv_mat=dummy_class,
                                 sigma = sigma,
                                 Int = TRUE)
}else{
  # Get average expression profile
  avg_expr <- CreateLibraryMatrix(t(dataDE@assays$RNA$data),
                                  dataDE@meta.data[,c("cell_ID","class")])
  # NicheDE object
  NDE_obj <- CreateNicheDEObject(counts_mat = t(dataDE@assays$RNA$data),
                                 coordinate_mat = dataDE@meta.data[,c("x","y")],
                                 library_mat = avg_expr[colnames(dummy_class),],
                                 deconv_mat=dummy_class,
                                 sigma = sigma,
                                 Int = FALSE)
}


# Calculate effective niche (Gaussian kernel)
NDE_obj <- CalculateEffectiveNiche(NDE_obj)
NDE_obj@effective_niche

#Perform Niche-DE
NDE_obj <- benchmark_and_return(niche_DE(NDE_obj,
                                         num_cores = ncores,outfile = "",C = 0, M = 0,
                                         gamma = 0, print = T, Int = (dist == "negbin"),
                                         batch = F,self_EN = T,G = 1))
save(NDE_obj, file = paste0(data_dir,"ResultsDistance/NicheDE/nicheDE_Microglia_",sigma,"_",dist,".Rdata"))
