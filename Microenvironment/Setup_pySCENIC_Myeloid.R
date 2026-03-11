library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)
library(SeuratObject)
library(SeuratDisk)
library(RcisTarget)

setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Macrophages/SCENIC")

seurat.rna <- readRDS('../Reanalysis_scRNA_macrophage_filtered.rds')
DimPlot(seurat.rna, raster = FALSE, label = T)
Idents(seurat.rna)
DefaultAssay(seurat.rna) <- "RNA"

dgem <- seurat.rna@assays$RNA@counts
dim(dgem)
head(colnames(dgem)) # columns => cell IDs

# Extract cell-level metadata
cell.info <- seurat.rna@meta.data

#############################
# Compute logical vector for every gene reflecting whether
# it has more than zero count per cell & and is detected in at 1 per thousand of cells in the global dataset
c <- floor(dim(dgem)[2]*.001); c  
nonzero <- dgem > 0L
keep_genes <- rowSums(as.matrix(nonzero)) >= c
filtered_dgem <- dgem[keep_genes, ]

# Further exclude mitochondrial genes
idx <- grep("^MT-", rownames(filtered_dgem))
rownames(filtered_dgem)[idx]
filtered_dgem <- filtered_dgem[-idx, ]
dim(filtered_dgem) 

#############################

write.table(t(as.matrix(dgem)), 
            'GPaggr_counts_mtx.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

####Using filtered genes from above
file.name <- "pySCENIC_RNA_Filtered.loom"
project.title <- "pHGG_Myeloid"
build_loom(
  file.name = file.name,
  dgem = filtered_dgem,
  title = project.title,
  genome = "human",
)

loom <- open_loom(file.path = file.name, mode = "r+")
finalize(loom = loom)
