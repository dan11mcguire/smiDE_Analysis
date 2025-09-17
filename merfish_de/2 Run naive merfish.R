#--------------------#
#- Run Naive Models -#
#--------------------#

data_dir <- "~MERFISH/"
library(tidyverse)
library(Seurat)
library(smiDE)
library(peakRAM)
library(future)
library(future.apply)
library(furrr)
library(DESeq2)
benchmark_and_return <- function(expr) {
  result <- NULL
  mem <- peakRAM(result <- eval(expr))
  list(result = result, memory = mem)
}

load( paste0(data_dir, "Allen/WholeBrain_seurat_distance.Rdata"))
load(paste0(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
colnames(dataDE@meta.data)[4] <- "cell_ID"
var <- "kernelMicroglia1000"
ncores <- 5

# smiDE -------------------------------------------------------------------
plan("multisession", workers = ncores)
dataDE@meta.data$DEvar <- dataDE@meta.data[,var]

## Naive ----
genes <- rownames(dataDE@assays$RNA$data)

#Negative Binomial
print("Negative Binomial Naive")

naivesmiDE_nb <- benchmark_and_return(quote(future_map(genes, function(gene_index) {
  smi_de(
    formula = ~ offset(log(nCount_RNA)) + DEvar,
    assay_matrix = dataDE@assays$RNA$counts,
    metadata = dataDE@meta.data,
    family = "nbinom2",
    targets = gene_index
  )
}, .options = furrr_options(seed = TRUE))
))
save(naivesmiDE_nb, file = paste0(data_dir,"ResultsDistance/NaivesmiDE/Naive_smiDE_nb_",var,".Rdata"))


# Gaussian
print("Gaussian Naive")
naivesmiDE_gaus <- benchmark_and_return(quote(future_map(genes, function(gene_index) {
  smi_de(
    formula = ~ DEvar,
    assay_matrix = dataDE@assays$RNA$data,
    metadata = dataDE@meta.data,
    family = "gaussian",
    targets = gene_index
  )
}, .options = furrr_options(seed = TRUE))
))
save(naivesmiDE_gaus, file = paste0(data_dir,"ResultsDistance/NaivesmiDE/Naive_smiDE_gaus_",var,".Rdata"))

# SemiNaive ----
options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
load(file.path(data_dir, "Allen/Cont_pre_de_WholeBrain_distance.Rdata"))
genes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1, ratio < 1) %>%pull(target)

# # Negative Binomial
print("Negative Binomial SemiNaive")
seminaivesmiDE_nb <- benchmark_and_return(quote(future_map(genes, function(gene_index) {
  smi_de(
    formula = ~ offset(log(nCount_RNA)) + DEvar + RankNorm(otherct_expr) ,
    assay_matrix = dataDE@assays$RNA$counts,
    metadata = dataDE@meta.data,
    family = "nbinom2",
    targets = gene_index,
    pre_de_obj = pre_de_allen$result,
    neighbor_expr_cell_type_metadata_colname = "class"
  )
}, .options = furrr_options(seed = TRUE))
))
save(seminaivesmiDE_nb, file = paste0(data_dir,"ResultsDistance/NaivesmiDE/semiNaive_smiDE_nb_",var,".Rdata"))

# Gaussian
print("Gaussian semiNaive")
seminaivesmiDE_gaus <- benchmark_and_return(quote(future_map(genes, function(gene_index) {
  smi_de(
    formula = ~ DEvar + RankNorm(otherct_expr),
    assay_matrix = dataDE@assays$RNA$data,
    metadata = dataDE@meta.data,
    family = "gaussian",
    targets = gene_index,
    pre_de_obj = pre_de_allen$result,
    neighbor_expr_cell_type_metadata_colname = "class"
  )
}, .options = furrr_options(seed = TRUE))
))
save(seminaivesmiDE_gaus, file = paste0(data_dir,"ResultsDistance/NaivesmiDE/semiNaive_smiDE_gaus_",var,".Rdata"))