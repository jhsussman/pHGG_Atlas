library(patchwork)
library(gridExtra)
library(Seurat)
library(data.table)

SplitByPatient = FALSE
nveg = 2000

## for tumor only
seurat.rna = readRDS('/mnt/isilon/tan_lab/oldridged/share/pHGG/SeuratObjects/snRNA_merged_SeuratObj_QCFiltered_DoubletFiltered_downsampled_rPCA_integrated_clustered_withInferCNVTumorNormal_withNeftelTypes.rds')
seurat.rna = subset(seurat.rna, site == 'CHOP')
seurat.rna$cellType[seurat.rna$cellType == 'Mature Neuron'] = 'Neuron'
seurat.rna$cellType[seurat.rna$cellType == 'Mature Oligodendrocyte'] = 'Oligodendrocyte'
seurat.rna$cellType[seurat.rna$cellType == 'Neuroglial'] = 'Neuroglia'
seurat.rna$cellType[seurat.rna$TumorNormal == 'Normal' & seurat.rna$cellType == 'Neuroglia'] = 'Neuroglia_Normal'

seurat.rna <- subset(seurat.rna, cellType == 'Neuroglia')

DefaultAssay(seurat.rna) <- 'RNA'
seurat.rna[['integrated']] <- NULL

seurat.rna$sampleID = sapply(seurat.rna$sampleID, function(x) gsub('7316-', '', x) )

if(SplitByPatient){
  seurat.list = SplitObject(seurat.rna, split.by = 'patient_id')
}else{
  seurat.list = SplitObject(seurat.rna, split.by = 'sampleID')
}
rm(seurat.rna)

seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
})

#seurat.list = lapply(seurat.list, FUN = SCTransform, method = 'glmGamPoi')
features <- SelectIntegrationFeatures(seurat.list, nfeatures = nveg)
#seurat.list = PrepSCTIntegration(seurat.list, anchor.features = features)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#seurat.list = lapply(seurat.list, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                         anchor.features = features, dims = 1:30, 
                                       reduction = "rpca", k.anchor = 5)
seurat.rna <- IntegrateData(anchorset = anchors,  dims = 1:30, k.weight = 50)
DefaultAssay(seurat.rna) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.rna <- ScaleData(seurat.rna, verbose = FALSE)

seurat.rna <- RunPCA(seurat.rna, verbose = FALSE)
seurat.rna <- RunUMAP(seurat.rna, reduction = "pca", dims = 1:30)

seuratPath = paste0('SeuratObjects/seurat_rna_Neuroglia_rpca_bySample_veg', nveg, '.rds')
plotPath = paste0('SeuratObjects/seurat_rna_Neuroglia_rpca_bySample_veg', nveg, '.pdf')

if(SplitByPatient){
  seuratPath = paste0('SeuratObjects/seurat_rna_Neuroglia_rpca_byPatient_veg', nveg, '.rds')
  plotPath = paste0('SeuratObjects/seurat_rna_Neuroglia_rpca_byPatient_veg', nveg, '.pdf')
}
p1 <- DimPlot(seurat.rna, reduction = "umap", group.by = "timepoint")
p2 <- DimPlot(seurat.rna, reduction = "umap", group.by = "sampleID")

saveRDS(seurat.rna, file = seuratPath)

ggsave(plot = marrangeGrob(list(p1, p2), nrow=1, ncol=1), 
       width = 8, height = 6, 
       filename = plotPath, 
       device = 'pdf')


