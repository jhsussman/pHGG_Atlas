## Introduction ##
# This script extracts counts and cell type/state labels from a seurat
# object and writes to tsv (to be used as sc reference for CIBERSORTx)

# Load Packages
library(tidyverse)
library(Seurat)
library(data.table)

# Set paths to data and output dir
base.dir <- "/mnt/isilon/tan_lab/tumultyj/MetaPathwayPaper/Brain/Glioma/pHGG/"
data.dir <- paste0(base.dir,"snRNA/")
R.obj.dir <- paste0(base.dir,"R_Objects/")
figure.dir <- paste0(base.dir,"Figures/")
decon.dir <- paste0(base.dir,"Deconvolution/data/")

# Load seurat objects
tumor.seurat <- readRDS(paste0(R.obj.dir,"snRNA_seurat_neuroglia.rds"))
full.seurat <- readRDS(paste0(R.obj.dir,"snRNA_seurat_allCells.rds"))

# Get normal cells from full seurat
norm.seurat <- full.seurat[,full.seurat$merged_cellType!="Neuroglial"]
rm(full.seurat)

# Combine endothelial and pericytes to vascular cells
merged_cellType_condensed <- norm.seurat$merged_cellType
merged_cellType_condensed <- gsub("Endothelial Cell", "Vascular Cell", merged_cellType_condensed)
merged_cellType_condensed <- gsub("Pericyte", "Vascular Cell", merged_cellType_condensed)
norm.seurat$merged_cellType_cond <- merged_cellType_condensed

# Set idents to be used in sc reference
norm.seurat <- SetIdent(norm.seurat, value="merged_cellType_cond")
tumor.seurat <- SetIdent(tumor.seurat, value="cell_state1")

# Get counts (CIBERSORTx will normalize)
norm.counts.mtx <- norm.seurat@assays$RNA@counts
tumor.counts.mtx <- tumor.seurat@assays$RNA@counts

## Create combined counts matrix to deconvolute bulkRNA
## with all cell types/states
counts.mtx <- cbind(norm.counts.mtx,tumor.counts.mtx)
# Filter out genes not expressed in any cells
exp.genes <- rowSums(counts.mtx)!=0
counts.mtx <- counts.mtx[exp.genes,]
# Convert to data frame
counts.label.df <- as.data.frame(as.matrix(counts.mtx)) %>%
  rownames_to_column(var="Genes")
# Create row with cell type/state idents
row.1 <- cbind(data.frame("Genes"=c("GeneSymbol")),
               t(data.frame(norm.seurat@active.ident)),
               t(data.frame(tumor.seurat@active.ident)))
row.names(row.1) <- NULL
# Bind using data.table (fastest)
counts.label.df <- data.table::rbindlist(list(row.1,counts.label.df))
# Write table
write.table(counts.label.df, file=paste0(decon.dir,"snRNA_cellState1CellTypeCondLabeled_counts.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


## Repeat above with only tumor cell states to deconvolute
## tumor cell line bulkRNA
counts.mtx <- tumor.counts.mtx
# Filter out genes not expressed in any cells
exp.genes <- rowSums(counts.mtx)!=0
counts.mtx <- counts.mtx[exp.genes,]
# Convert to data frame
counts.label.df <- as.data.frame(as.matrix(counts.mtx)) %>%
  rownames_to_column(var="Genes")
row.1 <- cbind(data.frame("Genes"=c("GeneSymbol")),
               t(data.frame(tumor.seurat@active.ident)))
row.names(row.1) <- NULL
# Bind using data.table
counts.df.label <- data.table::rbindlist(list(row.1,counts.df.label))
# Write table
write.table(counts.df.label, file=paste0(decon.dir,"snRNA_tumor_cellState1Labeled_counts.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
