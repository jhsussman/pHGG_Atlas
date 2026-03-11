library(ggplot2)
library(Seurat)
library(Signac)
library(ggrastr)
library(dittoSeq)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(infercnv)
library(paletteer)
options(future.globals.maxSize = 80000 * 1024^2)
setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/Multiome")

mapped_data = readRDS("Seurat_Multiome_Annotated.rds")

readClusters = read.table(paste0("InferCNV/infercnv_subclusters.observation_groupings.txt"), sep = ' ', header = T)
freq_table <- table(readClusters$Annotation.Group)
filtered_table <- freq_table[freq_table > 0]
valid_categories <- names(filtered_table)
filtered_data <- readClusters[readClusters$Annotation.Group %in% valid_categories, ]
mapped_data = AddMetaData(mapped_data, filtered_data)
subset1 = subset(mapped_data, subset = Annotation.Group %in% names(table(filtered_data$Annotation.Group)))
Idents(subset1) = "Annotation.Group"
DimPlot(subset1, group.by = "Annotation.Group", reduction = "ref.umap", raster = F, label = T, label.box = T, repel = T) + coord_fixed()
DimPlot(subset1, group.by = "Annotation.Group", reduction = "ref.umap", 
        raster = F, split.by = "Annotation.Group", ncol = 4) + coord_fixed()

rev(sort(filtered_table))
table(subset1$orig.ident, subset1$Annotation.Group)

mapped_data$FromInferCNV = mapped_data$Annotation.Group
mapped_data$Annotation.Color = NULL
mapped_data$Annotation.Group
table(mapped_data$orig.ident, mapped_data$predicted.celltype.proj)
table(mapped_data$FromInferCNV, mapped_data$orig.ident)

mapped_data$FromInferCNV = as.character(mapped_data$FromInferCNV)
Idents(mapped_data) = "FromInferCNV"
mapped_data = RenameIdents(mapped_data, 
                           '1' = "Tumor", '2' = "Tumor", '3' = "Tumor", '4' = "Tumor", 
                           '5' = "Tumor", '6' = "Tumor", '7' = "Tumor", '8' = "Normal", 
                           '9' = "Tumor", '10' = "Tumor", '11' = "Tumor", '12' = "Tumor", 
                           '13' = "Tumor", '14' = "Tumor", '15' = "Tumor", '16' = "Tumor", 
                           '17' = "Tumor", '18' = "Tumor", '19' = "Tumor", '20' = "Tumor", 
                           '21' = "Tumor", '22' = "Tumor", '23' = "Tumor", '24' = "Tumor", 
                           '25' = "Tumor", '26' = "Tumor", '27' = "Tumor", '28' = "Tumor", 
                           '29' = "Tumor", '30' = "Tumor")
mapped_data$TumorNormal = Idents(mapped_data)
DimPlot(mapped_data, group.by = "TumorNormal", reduction = "ref.umap", raster = F) + coord_fixed()
DimPlot(mapped_data, group.by = "TumorNormal", reduction = "ref.umap", raster = F, split.by = "TumorNormal") + coord_fixed()
#saveRDS(mapped_data, "Seurat_Multiome_Annotated.rds")

tumor_only = subset(mapped_data, subset=TumorNormal=="Tumor")

#Project neoplastic cells 
ng_rna <- readRDS("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Analysis_Revisions/Neuroglia_RNA_with_UMAP_Model.rds")
ng_rna <- AddMetaData(ng_rna, metadata=readRDS("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Analysis_Revisions/Metadata_cell_state1.rds"))

tumor_colors <- paletteer_d("ggthemes::Classic_10", n = 10)
tumor_colors <- as.vector(tumor_colors)
names(tumor_colors) <- c("GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                         "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")

table(tumor_only$orig.ident)

###############################################
#Reintegrate all together 
table(refquery$sampleID)
table(projected_rna$orig.ident)

refquery@meta.data <- refquery@meta.data %>%
  mutate(sampleID = case_when(
    refquery$orig.ident == "pHGG_2594-2-R1" ~ "2594",
    refquery$orig.ident == "pHGG_5831-R2" ~ "5831",
    refquery$orig.ident == "pHGG_8595" ~ "8595",
    TRUE ~ sampleID  
  ))
table(refquery$sampleID)

nveg = 2000
seurat.rna = refquery
DefaultAssay(seurat.rna) <- 'RNA'
seurat.rna[['integrated']] <- NULL
seurat.rna[['prediction.score.celltype.fan']] <- NULL
seurat.rna[['prediction.score.celltype.couterier']] <- NULL
seurat.rna[['prediction.score.celltype.proj']] <- NULL
seurat.rna[['prediction.score.phgg_refmap']] <- NULL
seurat.rna = JoinLayers(seurat.rna)
seurat.rna@assays$RNA@layers$scale.data = NULL
seurat.rna@assays$RNA@layers$scale.data.2 = NULL

seurat.rna.fixed = CreateSeuratObject(counts = CreateAssayObject(counts = GetAssayData(seurat.rna, assay = "RNA", layer = "counts")))
seurat.rna.fixed = AddMetaData(seurat.rna.fixed, metadata = seurat.rna@meta.data)

seurat.list = SplitObject(seurat.rna.fixed, split.by = 'sampleID')

#Log Normalization
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst")
})
features <- SelectIntegrationFeatures(seurat.list, nfeatures = nveg)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#Continue
anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                  anchor.features = features, dims = 1:30, 
                                  reduction = "rpca", k.anchor = 5)
seurat.rna <- IntegrateData(anchorset = anchors,  dims = 1:30, k.weight = 50)
DefaultAssay(seurat.rna) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.rna <- ScaleData(seurat.rna, verbose = TRUE)
seurat.rna <- RunPCA(seurat.rna, verbose = TRUE)
seurat.rna <- RunUMAP(seurat.rna, reduction = "pca", dims = 1:30)

seurat.rna = FindNeighbors(seurat.rna, dims = 1:30)
seurat.rna = FindClusters(seurat.rna, resolution = c(0.3))
#saveRDS(seurat.rna, 'All_Tumor_Integrated_Seurat_RNA.rds')

seurat.rna = readRDS("All_Tumor_Integrated_Seurat_RNA.rds")
levels = c("GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
  "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")
seurat.rna$cell_state1 = factor(seurat.rna$cell_state1, levels = levels)
DimPlot(seurat.rna, reduction = "umap", group.by = "cell_state1", raster = F) + coord_fixed()
DimPlot(seurat.rna, reduction = "umap", group.by = "integrated_snn_res.0.3", raster = F, label = T, label.box = T) + coord_fixed()

Idents(seurat.rna) = "integrated_snn_res.0.3"
seurat.rna = RenameIdents(seurat.rna, 
                          "0" = "Transition 1",
                          "1" = "AC-like",
                          "2" = "OPC/NPC-like", 
                          "3" = "Transition 2", 
                          "4" = "Cycling", 
                          "5" = "GPC-like", 
                          "6" = "MES-like", 
                          "7" = "NEU-like", 
                          "8" = "HSP+ Cells", 
                          "9" = "Cycling", 
                          "10" = "OC-like", 
                          "11" = "AC-like", 
                          "12" = "Transition 1", 
                          "13" = "NEU-like")
seurat.rna$cell_state1_redo = Idents(seurat.rna)
#saveRDS(seurat.rna, 'All_Tumor_Integrated_Seurat2.rds')

seurat.rna$cell_state1_redo = factor(seurat.rna$cell_state1_redo, levels = levels)
DimPlot(seurat.rna, reduction = "umap", group.by = "cell_state1_redo", raster = F, cols = tumor_colors, label = T) + coord_fixed()
source("../Figures/Figure_Functions.R")
create_umap_plot(seurat.rna, metadata_category = "cell_state1_redo", 
                 color_mapping = tumor_colors, width = 15, height = 15,
                 ptsize = 3.5, alpha = 1, filename = "Figures/UMAP_All_Reintegrated_MultiomeHighlight.pdf")
seurat.rna$cell_state1_redo = as.character(seurat.rna$cell_state1_redo)
seurat.rna$cell_state1_redo_m <- ifelse(!is.na(seurat.rna$cell_state1), "Original", seurat.rna$cell_state1_redo)
seurat.rna$cell_state1_redo_m = factor(seurat.rna$cell_state1_redo_m, levels = c(levels, "Original"))
tumor_colors_add = tumor_colors
tumor_colors_add[["Original"]] = "#F9F9F9"

create_umap_plot(subset(seurat.rna, subset=cell_state1_redo_m=="Original"), metadata_category = "cell_state1_redo_m", 
                 color_mapping = tumor_colors_add, width = 13, height = 13,
                 ptsize = 6, alpha = 1, filename = "Figures/UMAP_RNA_Reintegrated_MultiomeHighlight_Background.pdf")
create_umap_plot(subset(seurat.rna, subset=cell_state1_redo_m=="Original", invert = T), metadata_category = "cell_state1_redo_m", 
                 color_mapping = tumor_colors_add, width = 13, height = 13,
                 ptsize = 6, alpha = 1, filename = "Figures/UMAP_RNA_Reintegrated_MultiomeHighlight.pdf")

multiome_only = subset(seurat.rna, subset=orig.ident %in% c("pHGG_2594-2-R1", "pHGG_5831-R2", "pHGG_8595"))
p0 = dittoBarPlot(multiome_only, var = "cell_state1_redo", group.by = "orig.ident", retain.factor.levels = T, color.panel = tumor_colors)
ggsave(p0, filename = "Figures/Multiome_3Mapped_Bars.pdf", width = 3.25, height = 4.75)

#saveRDS(multiome_only, "Multiome_Only_RNA_States_Reintegrated.rds")
new_metadata = multiome_only@meta.data
#saveRDS(new_metadata, "Tumor_Multiome_Metadata.rds")

rownames(new_metadata) = gsub("^TRUE_", "", rownames(new_metadata))
projected_individual = AddMetaData(projected_individual, metadata = multiome_only@meta.data)

table(projected_individual$cell_state1_redo)
sum(table(projected_individual$cell_state1_redo))
projected_individual[["SCT"]] = NULL
projected_individual[["prediction.score.phgg_refmap"]] = NULL
projected_individual[["prediction.score.celltype.proj"]] = NULL

projected_individual = JoinLayers(projected_individual)
saveRDS(projected_individual, "Tumor_Multiome_Only_Seurat_Annotated.rds")

