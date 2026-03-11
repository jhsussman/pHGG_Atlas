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
#Get endothelial cells from original, and save
original_seurat <- readRDS("Integrated_Seurat_Objects/Integrated_Pix4_rPCA_3_50000_CLR2_NoCC_DAPI.RDS")
original_seurat$annotations_lv1_full <- readRDS("Metadata/Coarse_Annotations_Metadata.RDS")
endothelial_seurat <- subset(original_seurat, subset=annotations_lv1_full=="Endothelial Cells")
endothelial_seurat$annotations_lv2 <- "Endothelial Cells"
saveRDS(endothelial_seurat, file = "Annotated_Endothelial_10162023.RDS")

#Import RDS Files 
neuroglia_seurat <- readRDS("Annotated_Seurat_Objects/Annotated_Neuroglia_10162023.RDS")
neuroglia_seurat[["CODEX"]] <- JoinLayers(neuroglia_seurat[["CODEX"]])
neuroglia_seurat[["sketch"]] <- NULL
myeloid_seurat <- readRDS("Annotated_Seurat_Objects/Annotated_Myeloid_10162023.RDS")
myeloid_seurat[["CODEX"]] <- JoinLayers(myeloid_seurat[["CODEX"]])
myeloid_seurat[["sketch"]] <- NULL
tcells_seurat <- readRDS("Annotated_Seurat_Objects/Annotated_TCells_10162023.RDS")
tcells_seurat[["CODEX"]] <- JoinLayers(tcells_seurat[["CODEX"]])
endothelial_seurat <- readRDS("Annotated_Seurat_Objects/Annotated_Endothelial_10162023.RDS")
endothelial_seurat[["CODEX"]] <- JoinLayers(endothelial_seurat[["CODEX"]])
DefaultAssay(endothelial_seurat) <- "CODEX"
endothelial_seurat[["sketch"]] <- NULL

#Merge all together 
all_merged <- merge(neuroglia_seurat, y = list(myeloid_seurat, tcells_seurat, endothelial_seurat))
DefaultAssay(all_merged) <- "CODEX"
table(all_merged$annotations_lv2)
saveRDS(all_merged, "Annotated_Seurat_Objects/Merged_Annotated_10162023.RDS")

################
#START HERE
###############
#Remove cells after revising the ROI for sample 161, and one with 0 
all_merged <- readRDS("Annotated_Seurat_Objects/Merged_Annotated_10162023.RDS")
Idents(all_merged) <- "annotations_lv2"
all_merged[["CODEX"]] <- JoinLayers(all_merged[["CODEX"]])
seurat161 <- readRDS("Individual_Seurat_Objects/CPTCA_7316-161_Pix4_revised.RDS")
cells161 <- colnames(subset(all_merged, subset=orig.ident=="7316-161_Pix4"))
cellsnew161 <- paste0("pHGG_161_",colnames(seurat161))
cellstoremove <- setdiff(cells161, cellsnew161)
cellstoremoveall <- c(cellstoremove, "pHGG_5335_0")
all_merged_filtered <- subset(all_merged, cells = cellstoremoveall, invert = TRUE)
table(all_merged_filtered$annotations_lv2)
all_merged_filtered <- subset(all_merged_filtered, subset = annotations_lv2 %in% c("Artifact", "RBCs"), invert = TRUE)
table(all_merged_filtered$annotations_lv2)

#Clean it up 
all_merged_filtered[["sketch"]] = NULL
all_merged_filtered$leverage.score = NULL
VariableFeatures(all_merged_filtered) = important_features
all_merged_filtered <- NormalizeData(object = all_merged_filtered, normalization.method = "CLR", margin = 2)

#Use sketch integration
VariableFeatures(all_merged_filtered) = important_features
all_merged_filtered2 <- split(all_merged_filtered, f = all_merged_filtered$orig.ident)
all_merged_filtered2 <- SketchData(object = all_merged_filtered2, ncells = sketch_count, verbose = TRUE,
                                   method = "LeverageScore", sketched.assay = "sketch", over.write = TRUE)
DefaultAssay(all_merged_filtered2) <- "sketch"
all_merged_filtered2 <- ScaleData(all_merged_filtered2, verbose = T)
all_merged_filtered2 <- RunPCA(all_merged_filtered2, features = important_features, verbose = T)
saveRDS(all_merged_filtered2, "Annotated_Seurat_Objects/Merged_Annotated_Filtered_1_10162023.RDS")

#Integrate along sketched assay
all_merged_filtered2 <- readRDS("Annotated_Seurat_Objects/Merged_Annotated_Filtered_1_10162023.RDS")
integrated_seurat <- IntegrateLayers(object = all_merged_filtered2, method = RPCAIntegration, features = important_features,
                                     verbose = T, orig = "pca", group.by = "orig.ident", dims = 1:15,
                                     new.reduction = "integrated.rpca", layers = Layers(all_merged_filtered2, search = "data"))
DefaultAssay(integrated_seurat) <- "sketch"
VariableFeatures(integrated_seurat) <- important_features 

integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:15, reduction = "integrated.rpca", 
                             reduction.name = "umap.rpca", return.model = TRUE)
VariableFeatures(integrated_seurat)
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "orig.ident", raster = FALSE, alpha = 0.1) + coord_fixed()

#Reproject the full dataset onto this integration 
VariableFeatures(integrated_seurat) <- important_features 
VariableFeatures(integrated_seurat)
integrated_seurat <- ProjectIntegration(object = integrated_seurat, sketched.assay = "sketch", 
                                        features = important_features,
                                        assay = "CODEX", reduction = "integrated.rpca")
options(future.globals.maxSize = 8000 * 1024^2)
integrated_seurat <- ProjectData(object = integrated_seurat, sketched.assay = "sketch", assay = "CODEX", 
                                 sketched.reduction = "integrated.rpca", 
                                 full.reduction = "integrated.rpca.full", umap.model = 'umap.rpca', dims = 1:15)
saveRDS(integrated_seurat, "Annotated_Seurat_Objects/Merged_Annotated_Filtered_2_10162023.RDS")

#Plotting with all cells 
DefaultAssay(integrated_seurat) <- "CODEX"
DimPlot(integrated_seurat, reduction = "full.umap.rpca", label = T, label.box = T, 
        group.by = "annotations_lv2", raster = FALSE, alpha = 0.05) + coord_fixed() 
all_features = rownames(integrated_seurat)
RidgePlot(integrated_seurat, group.by = "annotations_lv2", 
          features = paste0("codex_",all_features[1:12]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "annotations_lv2", 
          features = paste0("codex_",all_features[13:24]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "annotations_lv2", 
          features = paste0("codex_",all_features[25:36]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "annotations_lv2", 
          features = paste0("codex_",all_features[37:48]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "annotations_lv2", 
          features = paste0("codex_",all_features[49:length(all_features)]), slot = 'data')

