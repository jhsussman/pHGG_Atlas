## Introduction ##
# This script loads and filters bulk RNA data downloaded from
# Cavatica and the OpenPedCan project

# Load Packages
library(tidyverse)
library(Seurat)

# Set paths to data and output dir
base.dir <- "/mnt/isilon/tan_lab/tumultyj/MetaPathwayPaper/Brain/Glioma/pHGG/"
data.dir <- paste0(base.dir,"bulkRNA/OpenPedCan/")
R.obj.dir <- paste0(base.dir,"R_Objects/")
output.dir <- paste0(base.dir,"Deconvolution/data/")
path.strings <- strsplit(data.dir,"/")[[1]]
save.prefix <- paste0(path.strings[length(path.strings)-1],"_",path.strings[length(path.strings)],"_")

# Read in gene counts
# OpenPedCan TPM file
counts.df <- readRDS(paste0(data.dir,"gene-expression-rsem-tpm-collapsed.rds"))
counts.df <- rownames_to_column(counts.df, var="gene_names")

# Read in clinical data to filter
clinical.data <- read_delim(paste0(data.dir,"histologies_v15_7-13-24.tsv"), delim = "\t")
# Include all PBTA data
clinical.data.filt <- subset(clinical.data, subset = cohort=="PBTA")
# Only solid tissue RNA-seq samples
clinical.data.filt <- subset(clinical.data.filt, experimental_strategy=="RNA-Seq" & 
                               RNA_library!="exome_capture" &
                               composition=="Solid Tissue" &
                               sample_type=="Tumor")
# Clean up clinical data?
empty.cols <- sapply(colnames(clinical.data.filt), function(cn) {
                       return(all(is.na(clinical.data.filt[[cn]])))
                     })
clinical.data.filt <- clinical.data.filt[,!empty.cols]

# Subset counts data frame to only include filtered biospecimens
counts.filt.df <- counts.df[,c(1,match(clinical.data.filt$Kids_First_Biospecimen_ID,
                                       colnames(counts.df)))]
# Remove genes with 0 expression across all samples
keep <- rowSums(counts.filt.df[,2:dim(counts.filt.df)[2]])!=0
counts.filt.df <- counts.filt.df[keep,]

# Write to a tsv file for use with CIBERSORTx
write.table(counts.filt.df, file=paste0(output.dir,save.prefix,"allPBTA_TPM.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


## Make seurat object for survival analysis
counts.mtx <- counts.filt.df[,-1]
rownames(counts.mtx) <- counts.filt.df[,1]
bulk.seurat <- CreateSeuratObject(counts.mtx, assay = "bulkRNA", min.cells=0, min.features=0)
meta.data <- as.data.frame(clinical.data.filt[,c(1,3,4,6,8,9,10,12,15:42)])
rownames(meta.data) <- pull(clinical.data.filt,2)
bulk.seurat <- AddMetaData(bulk.seurat, meta.data)
saveRDS(bulk.seurat, file=paste0(R.obj.dir,"bulkRNA_allPBTA_TPMasCounts.rds"))

## Cell line kallisto counts file
data.dir <- paste0(base.dir,"bulkRNA/cellLine/")
cell.line.counts <- read_delim(paste0(data.dir,"Counts_Cell_Lines_Kallisto.txt"))
colnames(cell.line.counts)[1] <- "gene_names"
# Remove genes with 0 expression across all samples
keep <- rowSums(cell.line.counts[,2:dim(cell.line.counts)[2]])!=0
cell.line.counts <- cell.line.counts[keep,]
# Check for NA genes?
sum(is.na(cell.line.counts[,1]))
keep <- !is.na(cell.line.counts[,1])
cell.line.counts <- cell.line.counts[keep,]

# Write to a tsv file for use with CIBERSORTx
write.table(cell.line.counts, file=paste0(output.dir,"bulkRNA_cellLine_Kallisto_counts.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
