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
ob.list <- c()
for(name in sample_names){
  ob.list[[name]] <- readRDS(paste0("Individual_Seurat_Objects/CPTCA_7316-",name,"_",pix,".RDS"))
}
ob.list
combined <- merge(x = ob.list[[1]], y = unlist(ob.list[-1]), add.cell.ids = paste0("pHGG_", sample_names))
combined <- JoinLayers(combined)
print(important_features)
VariableFeatures(combined) = important_features
combined <- NormalizeData(object = combined, normalization.method = "CLR", margin = 2)
saveRDS(combined, paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_0_CLR2_NoCC_DAPI.RDS"))

#Use sketch integration
combined <- readRDS(paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_0_CLR2_NoCC_DAPI.RDS"))
VariableFeatures(combined) = important_features
VariableFeatures(combined)
combined <- split(combined, f = combined$orig.ident)
combined <- SketchData(object = combined, ncells = as.numeric(sketch_count), method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(combined) <- "sketch"
combined <- ScaleData(combined, verbose = T)
combined <- RunPCA(combined, features = important_features, verbose = T)
#saveRDS(combined, paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_1_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))

#Integrate along sketched assay
integrated_seurat <- IntegrateLayers(object = combined, method = RPCAIntegration, features = important_features,
                                     verbose = T, orig = "pca", group.by = "orig.ident", dims = 1:15,
                                     new.reduction = "integrated.rpca", layers = Layers(combined, search = "data"))
DefaultAssay(integrated_seurat) <- "sketch"
VariableFeatures(integrated_seurat) <- important_features 
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.rpca", dims = 1:15)
integrated_seurat <- FindClusters(integrated_seurat, algorithm = 2, 
                                  resolution = c(0.3,0.4,0.5)) 
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:15, reduction = "integrated.rpca", 
                             reduction.name = "umap.rpca", return.model = TRUE)
saveRDS(integrated_seurat, paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_2_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))


#Some plots 
#integrated_seurat <- readRDS(paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_2_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))
VariableFeatures(integrated_seurat)
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "orig.ident", raster = FALSE, alpha = 0.1) + coord_fixed()
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "sketch_snn_res.0.5", label = T, 
        label.box = T, raster = FALSE, alpha = 0.1) + coord_fixed()
VlnPlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
        features = paste0("sketch_",important_features[1:18]), slot = 'data', pt.size = 0)
VlnPlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
        features = paste0("sketch_",important_features[19:37]), slot = 'data', pt.size = 0)
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
          features = paste0("sketch_",important_features[1:16]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
          features = paste0("sketch_",important_features[17:32]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "sketch_snn_res.0.4", 
          features = paste0("sketch_",important_features[33:37]), slot = 'data')


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
pdf(paste0("Marker_Plots_Integrated_",pix,"_",as.character(sketch_count),"_CLR2_NoCC_DAPI.pdf"))
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
#integrated_seurat <- readRDS("Integrated_Seurat_Objects/Integrated_Pix4_rPCA_3b_50000.RDS")
options(future.globals.maxSize = 8000 * 1024^2)
integrated_seurat <- ProjectData(object = integrated_seurat, sketched.assay = "sketch", assay = "CODEX", 
                                 sketched.reduction = "integrated.rpca", 
                                 refdata = list(snn_res.0.5_full = "sketch_snn_res.0.5", 
                                                snn_res.0.4_full = "sketch_snn_res.0.4",
                                                snn_res.0.3_full = "sketch_snn_res.0.3"
                                 ), 
                                 full.reduction = "integrated.rpca.full", umap.model = 'umap.rpca', dims = 1:15)
saveRDS(integrated_seurat, paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_3_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))

######
#Generate Masks
######
annotations_col <- "snn_res.0.4_full"
cluster_factor <- as.numeric(as.factor(integrated_seurat$snn_res.0.4_full))
integrated_seurat$cluster_levels <- cluster_factor
conversion_table <- data.frame(
  Level = levels(factor(cluster_factor)),
  Numeric = as.numeric(levels(factor(integrated_seurat$snn_res.0.4_full)))
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
  write.csv(annotations_df, file = paste0("Annotation_Masks/",name,"_",annotations_col,"_",pix,"_NoCC_DAPI.csv"))
}

#Define first level annotations
table(integrated_seurat$orig.ident, integrated_seurat$snn_res.0.4_full)
DefaultAssay(integrated_seurat) <- "CODEX"
Idents(integrated_seurat) <- "snn_res.0.4_full"
integrated_seurat <- RenameIdents(integrated_seurat, 
                                  "0" = "Neuroglia", 
                                  "1" = "Neuroglia",
                                  "2" = "Neuroglia",
                                  "3" = "Myeloid Cells",
                                  "4" = "Endothelial Cells",
                                  "5" = "Neuroglia",
                                  "6" = "Neuroglia",
                                  "7" = "Myeloid Cells",
                                  "8" = "Myeloid Cells",
                                  "9" = "Neuroglia",
                                  "10" = "T Cells",
                                  "11" = "Myeloid Cells",
                                  "12" = "Neuroglia",
                                  "13" = "Neuroglia",
                                  "14" = "Myeloid Cells",
                                  "15" = "Neuroglia",  
                                  "16"= "T Cells",
                                  "17"="Neuroglia")
integrated_seurat$annotations_lv1_full <- Idents(integrated_seurat)
table(integrated_seurat$annotations_lv1_full)

DimPlot(integrated_seurat, reduction = "full.umap.rpca", 
        group.by = "annotations_lv1_full", raster = FALSE, alpha = 0.05) + coord_fixed()
DimPlot(integrated_seurat, reduction = "full.umap.rpca", label = T, label.box = T, 
        group.by = "snn_res.0.4_full", raster = FALSE, alpha = 0.1) + coord_fixed()
metadata = integrated_seurat@meta.data
saveRDS(metadata$annotations_lv1_full, "Metadata/Coarse_Annotations_Metadata.RDS")
 
#Read and load image
#integrated_seurat <- readRDS(paste0("Integrated_Seurat_Objects/Integrated_",pix,"_rPCA_3_",as.character(sketch_count),"_CLR2_ALL.RDS"))
