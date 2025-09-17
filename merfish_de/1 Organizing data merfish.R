# ----------------------------------------#
# - Mouse Brain - Allen Institute --------#
# ----------------------------------------#
#https://alleninstitute.github.io/abc_atlas_access/notebooks/zhuang_merfish_tutorial.html
#https://github.com/AllenInstitute/abc_atlas_access/blob/main/descriptions/MERFISH-C57BL6J-638850.md
#https://knowledge.brain-map.org/data/LVDBJAW8BI5YSS1QUBG
#https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/Allen-CCF-2020/20230630/
#https://www-nature-com.offcampus.lib.washington.edu/articles/s41586-023-06812-z

library(readr)
library(dplyr)
library(tidyverse)
library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(Seurat)
library(peakRAM)
library(smiDE)
library(nicheDE)

# Paths
data_dir <- "~MERFISH/Allen/"

# Load files
# coordinates
coords <- read_csv(file.path(data_dir, "ccf_coordinates.csv"),
                   col_types = cols(cell_label = col_character()))
parcelannot <- read_csv(file.path(data_dir, "parcellation_to_parcellation_term_membership.csv"))%>%
  select(parcellation_index,#parcellation_term_name,
         parcellation_term_acronym,
         parcellation_term_set_name) %>%
  pivot_wider(names_from = parcellation_term_set_name, values_from =parcellation_term_acronym )
coords <- coords %>% left_join(parcelannot)

# meta data
metadata <- read_csv(file.path(data_dir, "cell_metadata.csv"),
                     col_types = cols(cell_label = col_character()))
clusterann <- read_csv(file.path(data_dir, "cluster_annotation_term.csv")) %>%
  merge(.,read_csv(file.path(data_dir, "cluster_to_cluster_annotation_membership.csv")),
        by.x = "label", by.y = "cluster_annotation_term_label", all.y = T) %>%
  select(cluster_alias,name, cluster_annotation_term_set_name.x) %>%
  pivot_wider(names_from = cluster_annotation_term_set_name.x, values_from =name )
metadata <- metadata %>% left_join(clusterann) %>%
  left_join(coords, by = "cell_label", suffix = c("","_ccf"))

# genes
genes <- read_csv(file.path(data_dir, "gene.csv"))

# counts
sce <- readH5AD(file.path(data_dir, "C57BL6J-638850-raw.h5ad"), reader="R")
sce <- sce[, colSums(assay(sce, "X"))  > 0]

## Convert to Seurat -------------------------------------------------------
cell_ids <- metadata %>% filter(brain_section_label %in% c("C57BL6J-638850.51")) %>%
  pull(cell_label)

metadata_filter <- as.data.frame(metadata[metadata$cell_label %in% cell_ids,])
rownames(metadata_filter) <- as.character(metadata_filter$cell_label)
sce_filter <- sce[, colnames(sce) %in% cell_ids]

seurat_allen <- CreateSeuratObject(
  counts = as(assay(sce_filter, "X"), "dgCMatrix"),
  meta.data = metadata_filter,
  assay = "RNA"
)
seurat_allen <- subset(seurat_allen, subset = nCount_RNA > 20)
seurat_allen <- NormalizeData(seurat_allen, normalization.method = "RC", verbose = FALSE)
genes <- rownames(seurat_allen@assays$RNA)[!grepl("Blank",rownames(seurat_allen@assays$RNA))]
seurat_allen <- subset(seurat_allen, features = genes)
save(seurat_allen, file = file.path(data_dir, "WholeBrain_seurat.Rdata"))


## Contamination metrics ---------------------------------------------------
benchmark_and_return <- function(expr) {
  result <- NULL
  mem <- peakRAM(result <- eval(expr))
  list(result = result, memory = mem)
}

seurat_allen <- subset(seurat_allen, subset = division == "Isocortex")
colnames(seurat_allen@meta.data)[4] <- "cell_ID"
overlap_metrics_allen <- benchmark_and_return((
  smiDE::overlap_ratio_metric(assay_matrix = seurat_allen@assays$RNA$counts
                              ,metadata = seurat_allen@meta.data
                              ,cellid_col = "cell_ID"
                              ,cluster_col = "class"
                              ,sdimx_col = "x"
                              ,sdimy_col = "y"
                              ,radius = 0.05
  )
))

pre_de_allen <- benchmark_and_return((
  pre_de(metadata = seurat_allen@meta.data,
         ref_celltype = "01 IT-ET Glut",
         cell_type_metadata_colname = "class",
         cellid_colname = "cell_ID",
         mm_radius = 0.05 ,
         split_neighbors_by_colname = "brain_section_label",verbose = T,
         counts = seurat_allen@assays$RNA$counts,
         adjacencies_only = FALSE,
         sdimx_col = "x",
         sdimy_col = "y"
  )
))
save(overlap_metrics_allen,pre_de_allen, file = file.path(data_dir, "Cont_pre_de_WholeBrain_distance.Rdata"))

# Getting data frame ----
colnames(seurat_allen@meta.data)[4] <- "cell_ID"

# Filtering gene
filteredgenes <- overlap_metrics_allen$result %>% filter(class == "01 IT-ET Glut") %>%
  filter(avg_cluster > 0.1) %>%pull(target)
seurat_allen <- subset(seurat_allen, features = filteredgenes)
seurat_allen <- subset(seurat_allen, subset = division == "Isocortex")
seurat_allen <- subset(seurat_allen, subset =
                      (class == "01 IT-ET Glut" | subclass == "334 Microglia NN"))

# Calculate effective niche ----
library(nicheDE)
dummy_class <- model.matrix(~class-1, data = seurat_allen@meta.data)
colnames(dummy_class) <- gsub("class","",colnames(dummy_class))

# Get average expression profile
avg_expr <- CreateLibraryMatrix(t(seurat_allen@assays$RNA$counts),
                                seurat_allen@meta.data[,c("cell_ID","class")])
# NicheDE object
NDE_obj <- CreateNicheDEObject(counts_mat = t(seurat_allen@assays$RNA$counts),
                               coordinate_mat = seurat_allen@meta.data[,c("x","y")],
                               library_mat = avg_expr[colnames(dummy_class),],
                               deconv_mat=dummy_class,
                               sigma = c(1000),
                               Int = TRUE)
NDE_obj <- CalculateEffectiveNiche(NDE_obj)
effecniche1000 <- data.frame(kernelMicroglia1000 = NDE_obj@effective_niche$`1000`[,2])
seurat_allen@meta.data <- cbind(seurat_allen@meta.data,effecniche500,effecniche1000)

# Filtering to neurons
dataDE <- subset(seurat_allen, subset = class == "01 IT-ET Glut")
dataDEmetadata <- dataDE@meta.data
dataDEcounts <- dataDE@assays$RNA$counts
dataDEnorm <- dataDE@assays$RNA$data

save(dataDE,file = paste0(data_dir, "WholeBrain_seurat_distance.Rdata"))
save(dataDEmetadata,dataDEcounts,dataDEnorm,file= paste0(data_dir, "WholeBrain_df_distance.Rdata"))