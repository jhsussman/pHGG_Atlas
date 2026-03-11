### This script takes as input the Seurat Objects from each CODEX sample in the atlas and creates a combined CODEX Atlas Object
library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr) 
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
source("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CODEX_Functions.R")
setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CPTCA_Full")

####SET VARIABLES####
sample_names <- c("7622", "6477", "5928", "4740", "4337", "3058", "942", "371", "339", "161", "5335")
pix = "Pix4"
sketch_count = 50000
important_features  <- c("DAPI", "CD31", "CD44","IBA1","NFP",
                         "CD4","ATRX","APOE","CD56","CD163",
                         "GFAP","CD11b","S100B","CD206","OLIG1",
                         "CD133","SOX2","Vimentin","MPO","HLA-DR",
                         "P2RY12","CD8","NeuN","Collagen IV", "CD14",
                         "CX3CR1","MOG","SPP1","CD3e","CD68","Nestin",
                         "CD16","TMEM119","OLIG2","GLUT1") 
####################
#Import RDS Files 
combined <- readRDS(paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_3_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))
combined$annotations_lv1_full <- readRDS("Metadata/Coarse_Annotations_Metadata.RDS")
table(combined$annotations_lv1_full)
DefaultAssay(combined) <- "CODEX"
combined0 <- subset(combined, subset=annotations_lv1_full=="Neuroglia")
combined0[["sketch"]] = NULL
combined0[["integrated.rpca"]]=NULL
combined0[["full.umap.rpca"]]=NULL
combined0[["integrated.rpca.full"]]=NULL
combined0[["umap.rpca"]]=NULL
combined0$leverage.score = NULL

combined0 <- JoinLayers(combined0)
VariableFeatures(combined0) = important_features
combined0 <- NormalizeData(object = combined0, normalization.method = "CLR", margin = 2)

#Use sketch integration
VariableFeatures(combined0) = important_features
combined0 <- split(combined0, f = combined0$orig.ident)
combined0 <- SketchData(object = combined0, ncells = sketch_count, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(combined0) <- "sketch"
combined0 <- ScaleData(combined0, verbose = T)
combined0 <- RunPCA(combined0, features = important_features, verbose = T)

#Integrate along sketched assay
integrated_seurat <- IntegrateLayers(object = combined0, method = RPCAIntegration, features = important_features,
                                     verbose = T, orig = "pca", group.by = "orig.ident", dims = 1:15,
                                     new.reduction = "integrated.rpca", layers = Layers(combined0, search = "data"))
DefaultAssay(integrated_seurat) <- "sketch"
VariableFeatures(integrated_seurat) <- important_features 
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.rpca", dims = 1:15)
integrated_seurat <- FindClusters(integrated_seurat, algorithm = 2, 
                                  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6)) 
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:15, reduction = "integrated.rpca", 
                             reduction.name = "umap.rpca", return.model = TRUE)
saveRDS(integrated_seurat, paste0("Neuroglia_Seurat_Objects/Neuroglia_",pix,"_rPCA_1_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))

#Some plots 
#integrated_seurat <- readRDS(paste0("Neuroglia_Seurat_Objects/Neuroglia_",pix,"_rPCA_1_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))
VariableFeatures(integrated_seurat)
other_markers <- c("Ki67", "PCNA", "PAX5", "PD-L1", "CD79a", "SOX4", 
                   "FOXP3", "CD47", "PD-1", "H3K27M", "PDGFRA", "ASCL1", "HIF1A", 
                   "p53", "CD38", "FOSL1", "EGFR")
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "orig.ident", raster = FALSE, alpha = 0.1) + coord_fixed()
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "sketch_snn_res.0.5", label = T, 
        label.box = T, raster = FALSE, alpha = 0.1) + coord_fixed()
VlnPlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
        features = paste0("sketch_",important_features[1:18]), slot = 'data', pt.size = 0)
VlnPlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
        features = paste0("sketch_",important_features[19:37]), slot = 'data', pt.size = 0)
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.5", 
          features = paste0("sketch_",important_features[1:16]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.5", 
          features = paste0("sketch_",important_features[17:32]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.5", 
          features = paste0("sketch_",important_features[33:35]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.5", 
          features = paste0("sketch_",other_markers[1:9]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.5", 
          features = paste0("sketch_",other_markers[10:17]), slot = 'data')

FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "sketch_HIF1A", 
            raster = FALSE, max.cutoff = 'q99', min.cutoff = 'q1', 
            alpha = 0.5, slot = "scale.data") + coord_fixed()
FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "Blank3", 
            raster = FALSE, max.cutoff = 'q99', min.cutoff = 'q1', 
            alpha = 0.5, slot = "data") + coord_fixed()
FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "Size", 
            raster = FALSE, #max.cutoff = 'q99', min.cutoff = 'q1', 
            alpha = 0.5, slot = "data") + coord_fixed()
FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "RowSum", 
            raster = FALSE, #max.cutoff = 'q99', min.cutoff = 'q1',
            alpha = 0.5, slot = "data") + coord_fixed() 
VlnPlot(integrated_seurat, group.by = "sketch_snn_res.0.5", 
        features = c("RowSum", "Size"), pt.size = 0)

#Print pdf  
pdf(paste0("Marker_Plots_Neuroglia_",pix,"_",as.character(sketch_count),"_CLR2_NoCC_DAPI.pdf"))
for(marker in important_features){
  print(marker)
  try(p <- FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = paste0("sketch_",marker), slot = 'scale.data',
                       raster = TRUE, raster.dpi = c(3000, 3000), max.cutoff = 'q99', min.cutoff = 'q1', pt.size = 3, alpha = 0.5) + 
        coord_fixed())
  #scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
  try(print(p))
} 
dev.off()
dev.off()

#Reproject the full dataset onto this integration 
VariableFeatures(integrated_seurat) <- important_features 
VariableFeatures(integrated_seurat)
integrated_seurat <- ProjectIntegration(object = integrated_seurat, sketched.assay = "sketch", 
                                        features = important_features,
                                        assay = "CODEX", reduction = "integrated.rpca")

#Including projection onto UMAP
options(future.globals.maxSize = 8000 * 1024^2)
integrated_seurat <- ProjectData(object = integrated_seurat, sketched.assay = "sketch", assay = "CODEX", 
                                 sketched.reduction = "integrated.rpca", 
                                 refdata = list(snn_res.0.5_full = "sketch_snn_res.0.5", 
                                                snn_res.0.4_full = "sketch_snn_res.0.4",
                                                snn_res.0.3_full = "sketch_snn_res.0.3"
                                 ), 
                                 full.reduction = "integrated.rpca.full", umap.model = 'umap.rpca', dims = 1:15)
saveRDS(integrated_seurat, paste0("Neuroglia_Seurat_Objects/Neuroglia_",pix,"_rPCA_2_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))

######
#Generate Masks
######
annotations_col <- "snn_res.0.5_full"
cluster_factor <- as.numeric(as.factor(integrated_seurat$snn_res.0.5_full))
integrated_seurat$cluster_levels <- cluster_factor
conversion_table <- data.frame(
  Level = levels(factor(cluster_factor)),
  Numeric = as.numeric(levels(factor(integrated_seurat$snn_res.0.5_full)))
)
conversion_table

for(name in sample_names){
  print(name)
  print("Subseting")
  sample_subset <- subset(integrated_seurat, subset=orig.ident==paste0("7316-",name,"_",pix))
  annotations_df <- as.data.frame(sample_subset$cluster_levels)
  rownames(annotations_df) <- sub(".*_", "", colnames(sample_subset)) 
  annotations_df$CellID <- sub(".*_", "", colnames(sample_subset)) 
  colnames(annotations_df) <- c("Annotation", "CellID")
  print("Writing Annotations")
  write.csv(annotations_df, file = paste0("Annotation_Masks/Neuroglia_",name,"_",annotations_col,"_",pix,"_NoCC_DAPI.csv"))
}

#Define second level annotations
DefaultAssay(integrated_seurat) <- "CODEX"
table(integrated_seurat$orig.ident, integrated_seurat$snn_res.0.5_full)
Idents(integrated_seurat) <- "snn_res.0.5_full"
integrated_seurat <- RenameIdents(integrated_seurat, 
                                  "0" = "Intermediate Tumor Cells", 
                                  "1" = "Intermediate Tumor Cells",
                                  "2" = "Proneural Tumor Cells",
                                  "3" = "Intermediate Tumor Cells",
                                  "4" = "Oligodendrocytes/White Matter",
                                  "5" = "Mesenchymal-1 Tumor Cells",
                                  "6" = "RBCs",
                                  "7" = "Proneural Tumor Cells",
                                  "8" = "Artifact",
                                  "9" = "Mesenchymal-2 Tumor Cells",
                                  "10" = "RBCs",
                                  "11" = "Intermediate Tumor Cells",
                                  "12" = "Mature Neurons",
                                  "13" = "Artifact",
                                  "14" = "Intermediate Tumor Cells",
                                  "15" = "Artifact", 
                                  "16" = "Intermediate Tumor Cells",
                                  "17" = "Artifact", 
                                  "18" = "Artifact",
                                  "19" = "Oligodendrocytes/White Matter",
                                  "20" = "Artifact",
                                  "21" = "Intermediate Tumor Cells",
                                  "22" = "Artifact",
                                  "23" = "Artifact")
integrated_seurat$annotations_lv2 <- Idents(integrated_seurat)
table(integrated_seurat$orig.ident, integrated_seurat$annotations_lv2)
saveRDS(integrated_seurat, "Annotated_Seurat_Objects/Annotated_Neuroglia_10162023.RDS")

DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "annotations_lv2", raster = FALSE, alpha = 0.1) + coord_fixed()

#Plotting with all cells 
DefaultAssay(integrated_seurat) <- "CODdEX"
DimPlot(integrated_seurat, reduction = "full.umap.rpca", label = T, label.box = T, 
        group.by = "snn_res.0.5_full", raster = FALSE, alpha = 0.1) + coord_fixed()

#Read latest file 
integrated_seurat <- readRDS("Neuroglia_Seurat_Objects/Neuroglia_Pix4_rPCA_2_50000_CLR2_NoCC_DAPI.RDS")
