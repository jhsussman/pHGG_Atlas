library(Seurat)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(readxl)
library(dplyr)
library(tibble)
library(data.table)
library(dittoSeq)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(ggpubr)
library(UCell)
library(car)
setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Macrophages")
source("../../Figures/Figure_Functions.R")

###Load latest file
seurat.rna <- readRDS(file = "../Final_Cohort_All_snRNA-seq.RDS")

#Subset myeloid
macrophage.rna <- subset(x=seurat.rna, subset = merged_cellType == "Macrophage/Microglia")
DimPlot(macrophage.rna, group.by = "timepoint")

#Add ribosomal gene content 
C<-GetAssayData(object = macrophage.rna, slot = "counts")
rb.genes <- rownames(macrophage.rna)[grep("^RP[SL]",rownames(macrophage.rna))]
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
macrophage.rna <- AddMetaData(macrophage.rna, percent.ribo, col.name = "percent.ribo")

DefaultAssay(macrophage.rna) <- 'RNA'
macrophage.rna <- SCTransform(macrophage.rna, method = "glmGamPoi", 
                              vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"),
                              verbose = TRUE)
macrophage.rna <- RunPCA(object = macrophage.rna)
macrophage.rna = RunUMAP(macrophage.rna, dims = 1:30, reduction = 'pca', reduction.name = 'umap')
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
macrophage.rna <- CellCycleScoring(macrophage.rna, s.features = s.genes, g2m.features = g2m.genes)

DimPlot(macrophage.rna, group.by = "patient_id") + coord_fixed()

#rPCA Integration
seurat.list = SplitObject(macrophage.rna, split.by = 'patient_id')
seurat.list = lapply(seurat.list, FUN = SCTransform, method = 'glmGamPoi', 
                     vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
features <- SelectIntegrationFeatures(seurat.list, nfeatures = 3000)
seurat.list = PrepSCTIntegration(seurat.list, anchor.features = features)
seurat.list = lapply(seurat.list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca", k.anchor = 20) 
macrophage.rna <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
macrophage.rna <- RunPCA(macrophage.rna, verbose = TRUE, reduction.name = "rpca_integrated")
macrophage.rna <- RunUMAP(macrophage.rna, reduction = "rpca_integrated", 
                          dims = 1:30, reduction.name = "umap_rpca")
#Timepoint
Idents(macrophage.rna) <- "timepoint"
macrophage.rna <- RenameIdents(macrophage.rna, "Initial CNS Tumor" = "Initial CNS Tumor", 
                               "Recurrence" = "Timepoint_2",
                               "Progressive (Non-Autopsy)" = "Timepoint_2", 
                               "Progressive (Autopsy)" = "Timepoint_2")
macrophage.rna$timemerge <- Idents(macrophage.rna)

#Clusters 
DefaultAssay(macrophage.rna) <- "integrated"
DimPlot(macrophage.rna, group.by = "timepoint") + coord_fixed()

macrophage.rna = FindNeighbors(macrophage.rna, reduction = 'rpca_integrated', 
                               dims = 1:30, verbose = T)
macrophage.rna = FindClusters(macrophage.rna, verbose = T,
                              resolution = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

#Annotation 
#TAM Gene lists
DefaultAssay(macrophage.rna) <- "RNA"
module.genes <- read.xlsx('TAM_S5.xlsx')
macrophage.rna <- AddModuleScore(macrophage.rna, features = list(na.omit(module.genes$MG_Markers)), name = "Microglia")
macrophage.rna <- AddModuleScore(macrophage.rna, features = list(module.genes$Mac_Markers), name = "Macrophage")

p1 <- FeaturePlot(macrophage.rna, features = 'Microglia1', min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca",
                  raster = FALSE) + ggtitle('Microglia') + coord_fixed() 
p2 <- FeaturePlot(macrophage.rna, features = 'Macrophage1',min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca", 
                  raster = FALSE) + ggtitle('Macrophage') + coord_fixed() 
p1+p2
saveRDS(macrophage.rna, file = "Reanalysis_scRNA_macrophage3.rds")


#DEGs
DefaultAssay(macrophage.rna) <- "RNA" 
Idents(macrophage.rna) <- "integrated_snn_res.0.2"
markers_res.0.2 <- FindAllMarkers(macrophage.rna, only.pos = TRUE,
                                  min.pct = 0.10, min.diff.pct = 0.15, logfc.threshold = 0.15)
macrophage.rna <- ScaleData(macrophage.rna)
markers_res.0.3 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top
write.table(markers_res.0.2, file = "Reanalysis_Macrophage_Markers_res0.2.txt", sep = '\t')
DotPlot(macrophage.rna, features = unique(top$gene), cols="RdBu") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) 
DoHeatmap(macrophage.rna, features = top$gene) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

#Projection to Harmonized Atlas 
harmonized_gbm_core <- readRDS("../Harmonized_Atlas/Harmonized_Reintegrated_rPCA_Author.RDS")

anchors <- FindTransferAnchors(
  reference = harmonized_gbm_core,
  query = macrophage.rna, 
  reference.assay = "RNA", query.assay = "RNA",
  reference.reduction = "integrated.rpca",
  dims = 1:50
)
Query_2_map <- MapQuery(
  anchorset = anchors,
  query = macrophage.rna.f,
  reference = harmonized_gbm_core,
  refdata = list(
    celltype.l1 = "annotation_level_1",
    celltype.l2 = "annotation_level_2",
    celltype.l3 = "annotation_level_3",
    celltype.l4 = "annotation_level_4"
  ),
  reference.reduction = "integrated.rpca", 
  reduction.model = "umap.rpca"
)
saveRDS(Query_2_map@meta.data, "Renalysis_Projection_2_Metadata.rds")
DimPlot(Query_2_map, group.by = "predicted.celltype.l4", label = TRUE,
        label.size = 3, repel = TRUE, reduction = 'umap_rpca') + coord_fixed()
FeaturePlot(Query_2_map, reduction = 'umap_rpca', features = "predicted.celltype.l1.score", cols = c("lightgrey", "darkred"),lbel.size = 3.5)
table(Query_2_map$integrated_snn_res.0.3, Query_2_map$predicted.celltype.l1)
dittoBarPlot(Query_2_map, "predicted.celltype.l3", group.by = "integrated_snn_res.0.6")

#################################################
#Remove putative neoplastic and contaminating cells and reintegrate 
macrophage.rna.f <- subset(macrophage.rna, integrated_snn_res.0.1 %in% c("1", "5"), invert = TRUE)
table(macrophage.rna.f$integrated_snn_res.0.3)
macrophage.rna.f[["SCT"]] = NULL
macrophage.rna.f[["integrated"]] = NULL

#rPCA Integration
seurat.list = SplitObject(macrophage.rna.f, split.by = 'patient_id')
seurat.list = lapply(seurat.list, FUN = SCTransform, method = 'glmGamPoi', 
                     vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
features <- SelectIntegrationFeatures(seurat.list, nfeatures = 3000)
seurat.list = PrepSCTIntegration(seurat.list, anchor.features = features)
seurat.list = lapply(seurat.list, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30,
                                  reduction = "rpca", k.anchor = 20) 
macrophage.rna.f <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
macrophage.rna.f <- RunPCA(macrophage.rna.f, verbose = TRUE, reduction.name = "rpca_integrated")
macrophage.rna.f <- RunUMAP(macrophage.rna.f, reduction = "rpca_integrated", 
                            dims = 1:30, reduction.name = "umap_rpca")
#Timepoint
Idents(macrophage.rna.f) <- "timepoint"
macrophage.rna.f <- RenameIdents(macrophage.rna.f, "Initial CNS Tumor" = "Initial CNS Tumor", 
                                 "Recurrence" = "Timepoint_2",
                                 "Progressive (Non-Autopsy)" = "Timepoint_2", 
                                 "Progressive (Autopsy)" = "Timepoint_2")
macrophage.rna.f$timemerge <- Idents(macrophage.rna.f)

#Clusters 
DefaultAssay(macrophage.rna.f) <- "integrated"
DimPlot(macrophage.rna.f, group.by = "timepoint") + coord_fixed()

macrophage.rna.f = FindNeighbors(macrophage.rna.f, reduction = 'rpca_integrated', 
                                 dims = 1:30, verbose = T)
macrophage.rna.f = FindClusters(macrophage.rna.f, verbose = T,
                                resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))

#Annotation 
#TAM Gene lists
DefaultAssay(macrophage.rna.f) <- "RNA"
module.genes <- read.xlsx('TAM_S5.xlsx')
macrophage.rna.f <- AddModuleScore(macrophage.rna.f, features = list(na.omit(module.genes$MG_Markers)), name = "Microglia")
macrophage.rna.f <- AddModuleScore(macrophage.rna.f, features = list(module.genes$Mac_Markers), name = "Macrophage")

p1 <- FeaturePlot(macrophage.rna.f, features = 'Microglia1', min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca",
                  raster = FALSE) + ggtitle('Microglia') + coord_fixed() 
p2 <- FeaturePlot(macrophage.rna.f, features = 'Macrophage1',min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca", 
                  raster = FALSE) + ggtitle('Macrophage') + coord_fixed() 
p1+p2
saveRDS(macrophage.rna.f, file = "Reanalysis_scRNA_macrophage_filtered2.rds")

#DEGs
DefaultAssay(macrophage.rna.f) <- "RNA" 
Idents(macrophage.rna.f) <- "integrated_snn_res.0.6"
markers_res.0.6 <- FindAllMarkers(macrophage.rna.f, only.pos = TRUE,
                                  min.pct = 0.10, min.diff.pct = 0.10, logfc.threshold = 0.10)
macrophage.rna.f <- ScaleData(macrophage.rna.f)
markers_res.0.6 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top
write.table(markers_res.0.6, file = "Reanalysis_Macrophage_Markers_res0.6.txt", sep = '\t')
DotPlot(macrophage.rna.f, features = unique(top$gene), cols="RdBu") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) 

#Signatures
signatures <- read_excel("S5_Signatures.xlsx")
DefaultAssay(macrophage.rna.f) <- "RNA" 
siglist = list()
siglist$Angiogenesis = as.character(na.omit(signatures$Angiogenesis))
macrophage.rna.f <- AddModuleScore(macrophage.rna.f, features = list(as.character(na.omit(signatures$M1))), name = "M1_Phenotype")
macrophage.rna.f <- AddModuleScore(macrophage.rna.f, features = list(as.character(na.omit(signatures$M2))), name = "M2_Phenotype")
macrophage.rna.f <- AddModuleScore(macrophage.rna.f, features = list(as.character(na.omit(signatures$Angiogenesis))), name = "Angiogenesis")
macrophage.rna.f <- AddModuleScore(macrophage.rna.f, features = list(as.character(na.omit(signatures$Phagocytosis))), name = "Phagocytosis")

p1 <- FeaturePlot(macrophage.rna.f, features = macrophage.rna.f, min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca",
                  raster = FALSE) + ggtitle('Angiogenesis') + coord_fixed() 
p2 <- FeaturePlot(macrophage.rna.f, features = 'VEGFA', min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca",
                  raster = FALSE) + ggtitle('Phagocytosis1') + coord_fixed() 
p3 <- FeaturePlot(macrophage.rna.f, features = 'M1_Phenotype1', min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca",
                  raster = FALSE) + ggtitle('M1_Phenotype') + coord_fixed() 
p4 <- FeaturePlot(macrophage.rna.f, features = 'M2_Phenotype1', min.cutoff = 'q1', max.cutoff = 'q90', reduction = "umap_rpca",
                  raster = FALSE) + ggtitle('M2_Phenotype') + coord_fixed() 
p1+p2+p3+p4
VlnPlot(macrophage.rna.f, features = "Angiogenesis1", group.by = 'integrated_snn_res.0.6', pt.size = 0)
VlnPlot(macrophage.rna.f, features = "Phagocytosis1", group.by = 'integrated_snn_res.0.6', pt.size = 0)

#Rename Clusters
Idents(macrophage.rna.f) <- "integrated_snn_res.0.6"
macrophage.rna.f <- RenameIdents(macrophage.rna.f, 
                                 `0`="Undetermined MG", 
                                 `1`="Pre-Active MG",
                                 `2`="BMD TAM 1",
                                 `3`="Lipid-Associated TAMs",
                                 `4`="Homeostatic MG",
                                 `5`="BMD TAM 2",
                                 `6`="BMD TAM 2",
                                 `7`="Proliferating Myeloid",
                                 `8`="Pro-Angiogenic TAM",
                                 `9`="BMD TAM 1",
                                 `10`="IFN-Responsive TAM",
                                 `11`="Inflammatory TAM",
                                 `12`="Pro-Angiogenic TAM",
                                 `13`="Dendritic Cells")
macrophage.rna.f$cell_labels1 <- Idents(macrophage.rna.f)
DimPlot(macrophage.rna.f, group.by = "cell_labels1", reduction = 'umap_rpca', 
        label = T) + coord_fixed() 
DimPlot(macrophage.rna.f, group.by = "timepoint", reduction = 'umap_rpca', 
        label = F) + coord_fixed() 
DimPlot(macrophage.rna.f, group.by = "patient_id", reduction = 'umap_rpca', 
        label = F) + coord_fixed() 
labels = c("Undetermined MG", "Pre-Active MG","BMD TAM 1","BMD TAM 2", "Lipid-Associated TAMs",
           "Homeostatic MG","Proliferating Myeloid","Pro-Angiogenic TAM","IFN-Responsive MG","Inflammatory TAM",
           "Dendritic Cells")
ggplot(macrophage.rna.f@meta.data, aes(x=cell_labels1, fill=patient_id)) + geom_bar(position = "fill") + scale_y_continuous(expand = c(0,0)) +
  theme_bw() + RotatedAxis() 

markers_celltypes1 <- FindAllMarkers(macrophage.rna.f, only.pos = TRUE,
                                  min.pct = 0.10, min.diff.pct = 0.10, logfc.threshold = 0.10)
markers_celltypes1 %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) -> top
write.table(markers_celltypes1, file = "Reanalysis_Macrophage_Markers_celltypes1.txt", sep = '\t')
DotPlot(macrophage.rna.f, features = unique(top$gene), cols="RdBu") + 
  RotatedAxis() + coord_flip() + 
  theme(axis.text.y = element_text(size = 6.5)) 

