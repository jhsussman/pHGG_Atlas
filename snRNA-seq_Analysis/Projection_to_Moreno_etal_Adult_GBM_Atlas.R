library(Seurat)
library(Matrix)
library(scales)
library(ggplot2)
library(stringr)
library(SingleCellExperiment)
library(anndata)
library(SeuratData)
library(SeuratDisk)
options(Seurat.object.assay.version = "v5")

setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Harmonized_Atlas")

harmonized_gbm_core <- readRDS(file = "Harmonized_GBM_Core.rds")
var.features.bool <- harmonized_gbm_core@assays$RNA@meta.features$highly.variable
features <-  as.character(harmonized_gbm_core@assays$RNA@meta.features$feature_name)
var.features <- c()
for(i in 1:length(var.features.bool)){
  if(var.features.bool[i]){
    a = as.character(features[i])
    var.features <- append(var.features, a)
  }
}
saveRDS(var.features, "Variable_Features.RDS")

counts <- GetAssayData(harmonized_gbm_core, assay = "RNA",slot = "counts") 
rownames(counts) <- harmonized_gbm_core@assays$RNA@meta.features$feature_name 
harmonized_new <- CreateSeuratObject(counts = counts, meta.data = harmonized_gbm_core@meta.data)


#Re-integration along dataset 
DefaultAssay(harmonized_new) <- "RNA"
harmonized_new <- NormalizeData(harmonized_new)
harmonized_new[["RNA"]] <- split(harmonized_new[["RNA"]], f = harmonized_new$author)
harmonized_new <- NormalizeData(harmonized_new)
harmonized_new <- FindVariableFeatures(harmonized_new, nfeatures = 3000)
harmonized_new <- ScaleData(harmonized_new)
harmonized_new <- RunPCA(harmonized_new)

harmonized_new <- IntegrateLayers(
  object = harmonized_new, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE
)
harmonized_new <- RunUMAP(harmonized_new, reduction = "integrated.rpca", 
                          dims = 1:30, reduction.name = "umap.rpca", return.model = TRUE)
saveRDS(harmonized_new, file = "Harmonized_Reintegrated_rPCA_Author.RDS")

#Read in harmonized atlas and project own data 
harmonized_new = readRDS(file = "Harmonized_Reintegrated_rPCA_Author.RDS")
phgg_snrna <- readRDS("../Final_Cohort_All_snRNA-seq.RDS")


anchors <- FindTransferAnchors(
  reference = harmonized_gbm_core,
  query = phgg_snrna, 
  features = var.features,
  reference.assay = "RNA", query.assay = "RNA",
  reference.reduction = "integrated.rpca",
  reduction = "pcaproject",
  dims = 1:50
)

phgg_mapped <- TransferData(anchorset = anchors, reference = harmonized_gbm_core, query = phgg_snrna,
                            refdata = list(celltype.l1 = "annotation_level_1", celltype.l2 = "annotation_level_2",
                                           celltype.l3 = "annotation_level_3", celltype.l4 = "annotation_level_4"))
saveRDS(phgg_mapped@meta.data, file = "Projected_Metadata.RDS")

phgg_mapped <- phgg_snrna
phgg_mapped@meta.data <- readRDS("Projected_Metadata_fixed.RDS")
DimPlot(phgg_mapped, group.by = "predicted.celltype.l1", label = TRUE,
        label.size = 3, raster = FALSE, alpha = 0.2) + ggtitle("Mapped Labels Level 1")
DimPlot(phgg_mapped, group.by = "predicted.celltype.l2", label = TRUE,
        label.size = 3, raster = FALSE, alpha = 0.2) + ggtitle("Mapped Labels Level 2")
DimPlot(phgg_mapped, group.by = "predicted.celltype.l3", label = TRUE,
        label.size = 3, raster = FALSE) + ggtitle("Mapped Labels Level 3")
DimPlot(phgg_mapped, group.by = "predicted.celltype.l4", label = TRUE,
        label.size = 3, raster = FALSE) + ggtitle("Mapped Labels Level 4")
