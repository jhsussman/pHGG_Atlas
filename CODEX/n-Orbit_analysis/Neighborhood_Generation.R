### This script takes as input the Seurat Objects from each CODEX sample in the atlas and finds neighborhoods for each sample individually
library(imcRtools)
library(Seurat)
library(tidyverse)
library(readr)
library(SingleCellExperiment)


options(future.globals.maxSize = 4e9)
options(Seurat.object.assay.version = "v5")
setwd("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma")


#Import RDS Files 
integrated_seurat <- readRDS("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CPTCA_Full/Annotated_Seurat_Objects/Merged_Annotated_Filtered_2_10162023.RDS")
integrated_seurat@meta.data%>%group_by(orig.ident) %>% summarize(count=n())


#Neighborhood analysis
samples = unique(integrated_seurat$orig.ident)
samples
for (sample in samples){
  sample_subset = subset(integrated_seurat, subset=orig.ident%in%c(sample))
  HGG.sce  <- SingleCellExperiment(sample_subset@assays$CODEX@layers$scale.data, colData=sample_subset@meta.data)
  names(colData(HGG.sce))[which(names(colData(HGG.sce))=="x.coord")]="Pos_X"
  names(colData(HGG.sce))[which(names(colData(HGG.sce))=="y.coord")]="Pos_Y"

  #Change the k to the number of nearest 
  HGG.sce <- buildSpatialGraph(HGG.sce, img_id = "orig.ident", type = "knn", k = 20)
  colPairNames(HGG.sce)
  HGG.sce <- aggregateNeighbors(HGG.sce, colPairName = "knn_interaction_graph", 
                              aggregate_by = "metadata", count_by = "annotations_lv2")

  #Change nstart to the number of clusters 
  cn_1 <- kmeans(HGG.sce$aggregatedNeighbors, centers = 15, nstart = 15, iter.max = 500)
  HGG.sce$cn_celltypes1 <- as.factor(cn_1$cluster)
  sample_subset$cn_celltypes1 <- as.factor(cn_1$cluster)
  df <- sample_subset@meta.data[, c("orig.ident", "x.coord", "y.coord", "annotations_lv2", "cn_celltypes1")]
  colnames(df) <- c("Sample_Name", "x.coord", "y.coord", "Cell_Type", "TCN_Label")
  write.csv(df, paste0("Nolan_", sample, ".csv"), row.names = FALSE)}