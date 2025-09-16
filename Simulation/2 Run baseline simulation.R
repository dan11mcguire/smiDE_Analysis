#----------------------------------------#
#-- Fit naive models to simulated data --#
#----------------------------------------#

dir.out = "~/Simulation/"

library(tidyverse)
library(Seurat)
library(smiDE)
library(peakRAM)
library(DESeq2)

benchmark_and_return <- function(expr) {
  expr <- substitute(expr)              
  calling_env <- parent.frame()         
  eval_env <- new.env(parent = calling_env)  
  mem <- peakRAM({
    eval_env$.__res__ <- eval(expr, envir = eval_env)
  })
  list(result = eval_env$.__res__, memory = mem)
}

# Read data ---------------------------------------------------------------

load(paste0(dir.out,paste0("SimulatedY.Rdata")))
n <- nrow(data)
data$cell_ID <- as.character(data$cell_ID)
rownames(data) <- data$cell_ID

# Negative Binomial ---------------------------------------------------------
p <- ncol(Ygaus)

## Seurat ------
# "negbinom", "poisson", "DESeq2"
seuratY <- CreateSeuratObject(counts = t(Ycounts),
                              data = t(Ygaus))
seuratY$Niche <- data$nichenum
Idents(seuratY) <- "Niche"
seuratY@meta.data$Niche <- as.factor(seuratY@meta.data$Niche)

seurat_negbinom <- benchmark_and_return((
  FindMarkers(
    seuratY,
    ident.1 = "1",
    ident.2 = "0",
    slot = "counts",
    fc.slot = "counts",
    test.use = "negbinom",
    min.cells.feature = -Inf,
    min.cells.group = -Inf,
    min.pct = -Inf, logfc.threshold = 0
  )
))


# Deseq2 
run_DEseq2 <- function(seuratY){
  dds <- DESeqDataSetFromMatrix(countData = seuratY@assays$RNA$counts,
                                colData = seuratY@meta.data, design = ~Niche)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  return(res)
}
seurat_deseq2 <- benchmark_and_return(run_DEseq2(seuratY))

# Gaussian ----------------------------------------------------------------

seurat_wilcox <- benchmark_and_return((
  FindMarkers(
    seuratY,
    ident.1 = "1",
    ident.2 = "0",
    slot = "data",
    test.use = "wilcox",
    min.cells.feature = -Inf,
    min.cells.group = -Inf,
    min.pct = -Inf, logfc.threshold = 0
  )
))

seurat_MAST <- benchmark_and_return((
  FindMarkers(
    seuratY,
    ident.1 = "1",
    ident.2 = "0",
    slot = "data",
    test.use = "MAST",
    min.cells.feature = -Inf,
    min.cells.group = -Inf,
    min.pct = -Inf, logfc.threshold = 0
  )
))

# Combining results
restodf <- function(seuratdata, modellabel, distlabel){
  df <- seuratdata$result %>%
    select(avg_log2FC, p_val, p_val_adj) %>%
    mutate(Time = seuratdata$memory$Elapsed_Time_sec,
           TotalRAM = seuratdata$memory$Total_RAM_Used_MiB,
           PeakRAM = seuratdata$memory$Peak_RAM_Used_MiB,
           model = modellabel,
           dist = distlabel
    )
  df$gene_id <- gsub("-","_",rownames(df))
  df
}

df <- seurat_deseq2$result %>% as.data.frame %>%
  dplyr::select(log2FoldChange, pvalue, padj)%>%
  dplyr::rename(avg_log2FC=log2FoldChange, p_val = pvalue, p_val_adj = padj) %>%
  mutate(Time = seurat_deseq2$memory$Elapsed_Time_sec,
         TotalRAM = seurat_deseq2$memory$Total_RAM_Used_MiB,
         PeakRAM = seurat_deseq2$memory$Peak_RAM_Used_MiB,
         model = "deseq2",
         dist = "negbin"
  ) %>% dplyr::add_rownames(var = "gene_id") %>%
  mutate(gene_id = gsub("-","_",gene_id)) %>%
  relocate(gene_id,.after=last_col()) %>%
  rbind(restodf(seurat_negbinom, "negbinom", "negbin")) %>%
  rbind(restodf(seurat_MAST, "MAST", "gaussian")) %>%
  rbind(restodf(seurat_wilcox, "wilcox", "gaussian")) %>%
  as.data.frame

results_naive <- df

save(results_naive,
     file = paste0(dir.out, "Results/PartialSimulation/ResultsNaiveModels.Rdata"))

