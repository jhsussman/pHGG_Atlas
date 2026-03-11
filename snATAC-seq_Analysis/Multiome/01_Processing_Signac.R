library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DoubletFinder)

setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/Multiome")

options(future.globals.maxSize = 80000 * 1024^2)
options(Seurat.object.assay.version = "v5")

#Load data
counts01 <- Read10X_h5("pHGG_2594-2-R1/outs/filtered_feature_bc_matrix.h5")
fragpath01 <- "pHGG_2594-2-R1/outs/atac_fragments.tsv.gz"
counts02 <- Read10X_h5("pHGG_5831-R2/outs/filtered_feature_bc_matrix.h5")
fragpath02 <- "pHGG_5831-R2/outs/atac_fragments.tsv.gz"
counts03 <- Read10X_h5("pHGG_8595/outs/filtered_feature_bc_matrix.h5")
fragpath03 <- "pHGG_8595/outs/atac_fragments.tsv.gz"

#Get annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

#Create a Seurat object containing the RNA adata
seurat01 <- CreateSeuratObject(
  counts = counts01$`Gene Expression`,
  assay = "RNA", project = "pHGG_2594-2-R1"
)
seurat02 <- CreateSeuratObject(
  counts = counts02$`Gene Expression`,
  assay = "RNA", project = "pHGG_5831-R2"
)
seurat03 <- CreateSeuratObject(
  counts = counts03$`Gene Expression`,
  assay = "RNA", project = "pHGG_8595"
)

#Create ATAC assay and add it to the object
seurat01[["ATACTemp"]] <- CreateChromatinAssay(
  counts = counts01$Peaks,
  sep = c(":", "-"),
  fragments = fragpath01,
  annotation = annotation
)
seurat02[["ATACTemp"]] <- CreateChromatinAssay(
  counts = counts02$Peaks,
  sep = c(":", "-"),
  fragments = fragpath02,
  annotation = annotation
)
seurat03[["ATACTemp"]] <- CreateChromatinAssay(
  counts = counts03$Peaks,
  sep = c(":", "-"),
  fragments = fragpath03,
  annotation = annotation
)

#QC check on ATAC, and RNA and filter 
DefaultAssay(seurat01) <- "ATAC"
seurat01 <- NucleosomeSignal(seurat01)
seurat01 <- TSSEnrichment(seurat01)
DefaultAssay(seurat01) <- "RNA"
seurat01[["percent.mt"]] <- PercentageFeatureSet(seurat01, pattern = "^MT-")
VlnPlot(
  object = seurat01,
  features = c("nCount_RNA", "nCount_ATACTemp", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0
)
seurat01_filtered <- subset(
  x = seurat01,
  subset = nCount_ATACTemp < 100000 &
    nCount_ATACTemp > 2000 &
    nFeature_RNA > 500 &
    nFeature_RNA < 10000 &
    TSS.enrichment > 2 & 
    percent.mt < 10
)

DefaultAssay(seurat02) <- "ATAC"
seurat02 <- NucleosomeSignal(seurat02)
seurat02 <- TSSEnrichment(seurat02)
DefaultAssay(seurat02) <- "RNA"
seurat02[["percent.mt"]] <- PercentageFeatureSet(seurat02, pattern = "^MT-")
VlnPlot(
  object = seurat02,
  features = c("nCount_RNA", "nCount_ATACTemp", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0
)
seurat02_filtered <- subset(
  x = seurat02,
  subset = nCount_ATACTemp < 100000 &
    nCount_ATACTemp > 2000 &
    nFeature_RNA > 500 &
    nFeature_RNA < 10000 &
    TSS.enrichment > 2 & 
    percent.mt < 10
)

DefaultAssay(seurat03) <- "ATAC"
seurat03 <- NucleosomeSignal(seurat03)
seurat03 <- TSSEnrichment(seurat03)
DefaultAssay(seurat03) <- "RNA"
seurat03[["percent.mt"]] <- PercentageFeatureSet(seurat03, pattern = "^MT-")
VlnPlot(
  object = seurat03,
  features = c("nCount_RNA", "nCount_ATACTemp", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
  ncol = 4,
  pt.size = 0
)
seurat03_filtered <- subset(
  x = seurat03,
  subset = nCount_ATACTemp < 100000 &
    nCount_ATACTemp > 2000 &
    nFeature_RNA > 500 &
    nFeature_RNA < 10000 &
    TSS.enrichment > 2 & 
    percent.mt < 10
)

FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.075, sct = FALSE){
  sweep.res.list <- paramSweep(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25,
                              pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                              sct = sct)
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25, 
                              pK = 0.09, nExp = nExp_poi.adj,
                              reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                              sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  return(seurat.rna)
}

seurat_all = list(seurat01, seurat02, seurat03)
saveRDS(seurat_all, "Seurat_All_List.rds")
seurat_filtered = list(seurat01_filtered, seurat02_filtered, seurat03_filtered)
saveRDS(seurat_filtered, "Seurat_Filtered_List.rds")

for (i in 1:length(seurat_filtered)) {
  DefaultAssay(seurat_filtered[[i]]) <- "RNA"
  seurat_filtered[[i]] <- SCTransform(seurat_filtered[[i]])
  seurat_filtered[[i]] <- RunPCA(seurat_filtered[[i]], features = VariableFeatures(object = seurat_filtered[[i]]))
  seurat_filtered[[i]] <- FindNeighbors(seurat_filtered[[i]], dims = 1:30)
  seurat_filtered[[i]] <- FindClusters(seurat_filtered[[i]], resolution = 1)
  seurat_filtered[[i]] <- FindDoublets(seurat_filtered[[i]], PCs = 1:30, sct = TRUE, exp_rate = 0.075)
}
saveRDS(seurat_filtered, "Seurat_Filtered_List_DoubletsMarked.rds")


#####################
#Recall features 
####################
seurat_filtered = readRDS("Seurat_Filtered_List_DoubletsMarked.rds")
signac.obj <- readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/SeuratObjects/seurat_atac_tumor_withTFMotifs.rds') ## with fragments and motif inf

granges_object = signac.obj@assays$ATAC@ranges

frags01 <- seurat_filtered[[1]]@assays$ATACTemp@fragments
frags02 <- seurat_filtered[[2]]@assays$ATACTemp@fragments
frags03 <- seurat_filtered[[3]]@assays$ATACTemp@fragments

seurat01_recounts <- FeatureMatrix(
  fragments = frags01,
  features = granges_object, process_n = 50000, 
  cells = colnames(seurat_filtered[[1]])
)
seurat02_recounts <- FeatureMatrix(
  fragments = frags02,
  features = granges_object, process_n = 50000, 
  cells = colnames(seurat_filtered[[2]])
)
seurat03_recounts <- FeatureMatrix(
  fragments = frags03,
  features = granges_object, process_n = 50000, 
  cells = colnames(seurat_filtered[[3]])
)

seurat01_filtered = seurat_filtered[[1]]
seurat02_filtered = seurat_filtered[[2]]
seurat03_filtered = seurat_filtered[[3]]
seurat01_filtered[["ATAC"]] <- CreateChromatinAssay(seurat01_recounts, fragments = frags01)
seurat02_filtered[["ATAC"]] <- CreateChromatinAssay(seurat02_recounts, fragments = frags02)
seurat03_filtered[["ATAC"]] <- CreateChromatinAssay(seurat03_recounts, fragments = frags03)

seurat01_filtered[["ATACTemp"]] = NULL 
seurat02_filtered[["ATACTemp"]] = NULL
seurat03_filtered[["ATACTemp"]] = NULL
seurat_filtered_merged = merge(seurat01_filtered, y = c(seurat02_filtered, seurat03_filtered))
saveRDS(seurat_filtered_merged, "Seurat_Multiome_All_DoubletsMarked.rds")
saveRDS(seurat_filtered_merged@meta.data, "Seurat_Multiome_All_Metadata.rds")

seurat_filtered_merged_singlet = subset(seurat_filtered_merged, subset=Doublet_Singlet=="Singlet")
saveRDS(seurat_filtered_merged_singlet, "Seurat_Multiome_All_DoubletsRemoved.rds")


