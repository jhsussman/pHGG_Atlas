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
sketch_count="NoSketch"
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
combined0 <- subset(combined, subset=annotations_lv1_full=="T Cells")
combined0[["sketch"]] = NULL
combined0[["integrated.rpca"]]=NULL
combined0[["full.umap.rpca"]]=NULL
combined0[["integrated.rpca.full"]]=NULL
combined0[["umap.rpca"]]=NULL
combined0$leverage.score = NULL 

combined0 <- JoinLayers(combined0)
VariableFeatures(combined0) = important_features
combined0 <- NormalizeData(object = combined0, normalization.method = "CLR", margin = 2)

VariableFeatures(combined0) = important_features
combined0 <- split(combined0, f = combined0$orig.ident)
combined0 <- ScaleData(combined0, verbose = T)
combined0 <- RunPCA(combined0, features = important_features, verbose = T)

#Integrate along all
integrated_seurat <- IntegrateLayers(object = combined0, method = RPCAIntegration, features = important_features,
                                     verbose = T, orig = "pca", group.by = "orig.ident", dims = 1:15,
                                     new.reduction = "integrated.rpca", layers = Layers(combined0, search = "data"))
VariableFeatures(integrated_seurat) <- important_features 
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.rpca", dims = 1:15)
integrated_seurat <- FindClusters(integrated_seurat, algorithm = 2, 
                                  resolution = c(0.01,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.6)) 
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:15, reduction = "integrated.rpca", 
                             reduction.name = "umap.rpca", return.model = TRUE)
saveRDS(integrated_seurat, paste0("TCells_Seurat_Objects/TCells_",pix,"_rPCA_1_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))

#Some plots 
#integrated_seurat <- readRDS(paste0("TCells_Seurat_Objects/TCells_",pix,"_rPCA_1_",as.character(sketch_count),"_CLR2_NoCC_DAPI.RDS"))
table(integrated_seurat$orig.ident, integrated_seurat$CODEX_snn_res.0.075)
other_markers <- c("Ki67", "PCNA", "PAX5", "PD-L1", "CD79a", "SOX4", 
                   "FOXP3", "CD47", "PD1", "H3K27M", "PDGFRA", "ASCL1", "HIF-1A", 
                   "p53", "CD38", "FOSL1", "EGFR", "Blank3")
VariableFeatures(integrated_seurat)
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "orig.ident", raster = FALSE) + coord_fixed()
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "CODEX_snn_res.0.075", label = T, 
        label.box = T, raster = FALSE, alpha = 0.8) + coord_fixed()
VlnPlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
        features = paste0("codex_",important_features[1:18]), slot = 'data', pt.size = 0)
VlnPlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
        features = paste0("codex_",important_features[19:37]), slot = 'data', pt.size = 0)
RidgePlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
          features = paste0("codex_",important_features[1:16]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
          features = paste0("codex_",important_features[17:32]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
          features = paste0("codex_",important_features[33:35]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
          features = paste0("codex_",important_features[33:35]), slot = 'data')
RidgePlot(integrated_seurat, group.by = "CODEX_snn_res.0.075", 
          features = paste0("codex_",important_features[33:35]), slot = 'data')

FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "codex_CD3e", 
            raster = FALSE, max.cutoff = 'q99', min.cutoff = 'q1', 
            alpha = 1, slot = "scale.data") + coord_fixed()
FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "Blank3", 
            raster = FALSE, max.cutoff = 'q99', min.cutoff = 'q1', 
            alpha = 0.5, slot = "data") + coord_fixed()
FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "Size", 
            raster = FALSE, #max.cutoff = 'q99', min.cutoff = 'q1', 
            alpha = 0.5, slot = "data") + coord_fixed()
FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = "RowSum", 
            raster = FALSE, #max.cutoff = 'q99', min.cutoff = 'q1',
            alpha = 0.5, slot = "data") + coord_fixed() 
VlnPlot(integrated_seurat, group.by = "CODEX_snn_res.0.5", 
        features = c("RowSum", "Size"), pt.size = 0)

#Print pdf  
pdf(paste0("Marker_Plots_TCells_",pix,"_",as.character(sketch_count),"_CLR2_NoCC_DAPI.pdf"))
for(marker in important_features){
  print(marker)
  try(p <- FeaturePlot(integrated_seurat, reduction = "umap.rpca", features = paste0("codex_",marker), slot = 'scale.data',
                       raster = TRUE, raster.dpi = c(3000, 3000), max.cutoff = 'q99', min.cutoff = 'q1', pt.size = 5, alpha = 1) + 
        coord_fixed())
  #scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu")))
  try(print(p))
} 
dev.off()
dev.off()

######
#Generate Masks
######
annotations_col <- "TSubtype"
cluster_factor <- as.factor(integrated_seurat$TSubtype)
integrated_seurat$cluster_levels <- cluster_factor
conversion_table <- data.frame(
  Level = levels(factor(cluster_factor)),
  Numeric = levels(factor(as.numeric(factor(integrated_seurat$TSubtype))))
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
  write.csv(annotations_df, file = paste0("Annotation_Masks/TCells_",name,"_",annotations_col,"_",pix,"_NoCC_DAPI.csv"))
}

###############
#Do annotation based on manual threshold 
#############

#Subset out artifact 
integrated_seurat = orig
combined0 <- subset(integrated_seurat, subset=CODEX_snn_res.0.075 %in% c("0", "1", "2"))
dim(combined0)
important_features  <- c("CD3e", "CD4", "CD8")
combined0 <- JoinLayers(combined0)
combined0 <- NormalizeData(object = combined0, normalization.method = "LogNormalize")
VariableFeatures(combined0) = important_features
RidgePlot(combined0, features = c("codex_CD3e", "codex_CD4", "codex_CD8"), 
          group.by = "annotations_lv1_full")
CD8Tcell_refined_ids <- colnames(subset(combined0, subset = codex_CD8 > 4.2))
combined0$TSubtype <- "CD4+ T Cell"
combined0$TSubtype[CD8Tcell_refined_ids] <- "CD8+ T Cell"
RidgePlot(combined0, features = c("codex_CD3e", "codex_CD4", "codex_CD8"), 
          group.by = "TSubtype")
integrated_seurat = combined0
saveRDS(integrated_seurat, paste0("TCells_Seurat_Objects/TCells_",pix,"_Manual.RDS"))


integrated_seurat$annotations_lv2 <- integrated_seurat$TSubtype
table(integrated_seurat$orig.ident, integrated_seurat$annotations_lv2)
DimPlot(integrated_seurat, reduction = "umap.rpca", 
        group.by = "annotations_lv2", raster = FALSE, alpha = 0.1) + coord_fixed()
saveRDS(integrated_seurat, "Annotated_Seurat_Objects/Annotated_TCells_10162023.RDS")

#Plotting with all cells 
DefaultAssay(integrated_seurat) <- "CODEX"
DimPlot(integrated_seurat, reduction = "full.umap.rpca", label = T, label.box = T, 
        group.by = "snn_res.0.4_full", raster = FALSE, alpha = 0.1) + coord_fixed()

#Read latest file 
integrated_seurat <- readRDS("TCells_Seurat_Objects/TCells_Pix4_Manual.RDS")

integrated_seurat <- readRDS("TCells_Seurat_Objects/TCells_Pix4_rPCA_1_50000_CLR2_NoCC_DAPI.RDS")
dim(integrated_seurat)
