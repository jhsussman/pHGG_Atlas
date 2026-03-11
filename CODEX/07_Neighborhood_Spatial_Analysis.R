### This script takes as input the Seurat Objects from each CODEX sample in the atlas and creates a combined CODEX Atlas Object
library(Seurat)
library(tidyverse)
library(patchwork)
library(readr)
library(dplyr) 
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(readr)
library(imcRtools)
library(SingleCellExperiment)
library(scales)
library(viridis)
library(ggpubr)
library(ggrastr)
library(pheatmap)
library(scrabbitr)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
setwd("/mnt/isilon/tan_lab/sussmanj/CODEX/Seurat_Analysis/CPTCA_Full")

####SET VARIABLES####
sample_names <- c("7622", "6477", "5928", "4740", "4337", "3058", "942", "371", "339", "161", "5335")
pix = "Pix4"
#Import RDS Files 
#Integrate along sketched assay
integrated_seurat <- readRDS("Annotated_Seurat_Objects/Merged_Annotated_Filtered_2_10162023.RDS")

######
#Generate Masks
######
annotations_col <- "annotations_lv2"
cluster_factor <- as.numeric(as.factor(integrated_seurat$annotations_lv2))
integrated_seurat$cluster_levels <- cluster_factor
conversion_table <- data.frame(
  Level = levels(as.factor(integrated_seurat$annotations_lv2)),
  Numeric = levels(factor(as.numeric(factor(cluster_factor))))
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
  write.csv(annotations_df, file = paste0("Annotation_Masks/Filtered_",name,"_",annotations_col,"_",pix,"_10172023.csv"))
}

#Heatmap of marker expression 
integrated_seurat_CODEX <- integrated_seurat
DefaultAssay(integrated_seurat_CODEX) <- "CODEX"
integrated_seurat_CODEX[["sketch"]] <- NULL

cluster.averages <- AverageExpression(integrated_seurat_CODEX, slot = "data", features = rownames(integrated_seurat), 
                                      return.seurat=FALSE, group.by = "annotations_lv2") 
cluster.averages <- as.data.frame(cluster.averages)
pheatmap(cluster.averages, cluster_cols = T, cluster_rows = T, scale = "row",
         color=colorRampPalette(c("dark blue", "white", "dark red"))(50), angle_col = 45)


#Neighborhood analysis using IMC Tools 
HGG.sce <- readRDS("Neighborhood_Files/HGG.sce.rds")
Idents(integrated_seurat) <- "annotations_lv2"
colors = DiscretePalette(16, palette = "alphabet", shuffle = FALSE)
table = integrated_seurat@meta.data
table$orig.ident <- factor(table$orig.ident)
rois = as.vector(unique(table$orig.ident))

ggplot(table, aes(x=orig.ident, fill=annotations_lv2), cols=colors) + geom_bar(position = "fill") + scale_y_continuous(expand = c(0,0)) +
  theme_bw() + RotatedAxis() + scale_fill_manual(values = colors)
theme(axis.text.x = element_text(angle = 90))


#Neighborhood analysis
HGG.sce  <- SingleCellExperiment(integrated_seurat@assays$CODEX@layers$scale.data, colData=integrated_seurat@meta.data)
names(colData(HGG.sce))[which(names(colData(HGG.sce))=="x.coord")]="Pos_X"
names(colData(HGG.sce))[which(names(colData(HGG.sce))=="y.coord")]="Pos_Y"

#Change the k to the number of nearest 
HGG.sce <- buildSpatialGraph(HGG.sce, img_id = "orig.ident", type = "knn", k = 20)
colPairNames(HGG.sce)
HGG.sce <- aggregateNeighbors(HGG.sce, colPairName = "knn_interaction_graph", 
                              aggregate_by = "metadata", count_by = "annotations_lv2")

#Change nstart to the number of clusters 
cn_1 <- kmeans(HGG.sce$aggregatedNeighbors, centers = 15, nstart = 50, iter.max = 500)
HGG.sce$cn_celltypes1 <- as.factor(cn_1$cluster)
integrated_seurat$cn_celltypes1 <- as.factor(cn_1$cluster)
for_plot <- prop.table(table(HGG.sce$cn_celltypes1, HGG.sce$annotations_lv2), 
                       margin = 1)
pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")
pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "row")
saveRDS(HGG.sce$cn_celltypes1, "Neighborhood_Files/CN_CellTypes1_k20_n15_backup.RDS")

##################
#Review and annotation 
sample_order <- c("7316-3058_Pix4","7316-339_Pix4","7316-942_Pix4",
                  "7316-161_Pix4","7316-5335_Pix4","7316-6477_Pix4",
                  "7316-371_Pix4", "7316-5928_Pix4","7316-7622_Pix4",
                  "7316-4337_Pix4","7316-4740_Pix4")
HGG.sce <- readRDS("Neighborhood_Files/HGG.sce.rds")
HGG.sce$cn_celltypes1 <- readRDS("Neighborhood_Files/CN_CellTypes1_k20_n15.RDS")
HGG.sce$cn_celltypes1 <- factor(HGG.sce$cn_celltypes1, levels = 1:15)
names(HGG.sce$cn_celltypes1) <- colnames(HGG.sce)
table(HGG.sce$cn_celltypes1)
integrated_seurat$cn_celltypes1 <- HGG.sce$cn_celltypes1[colnames(integrated_seurat)]
table(integrated_seurat$cn_celltypes1)

Idents(integrated_seurat) <- "cn_celltypes1"
integrated_seurat <- RenameIdents(integrated_seurat, `3`="CN1", 
                                  `2`="CN2",
                                  `5`="CN3",
                                  `13`="CN4",
                                  `14`="CN5",
                                  `4`="CN6",
                                  `11`="CN7",
                                  `8`="CN8",
                                  `1`="CN9",
                                  `10`="CN10",
                                  `15`="CN11",
                                  `7`="CN12",
                                  `9`="CN13", 
                                  `6`="CN14", 
                                  `12`="CN15")
integrated_seurat$neighborhood_numbers <- Idents(integrated_seurat)

Idents(integrated_seurat) <- "cn_celltypes1"
integrated_seurat <- RenameIdents(integrated_seurat, `3`="Immune enriched", 
                                  `2`="Gray matter",
                                  `5`="MES-2 enriched",
                                  `13`="Mixed/artifact",
                                  `14`="Perivascular",
                                  `4`="White matter",
                                  `11`="Infiltrating tumor",
                                  `8`="Proneural enriched",
                                  `1`="Tumor/macrophage",
                                  `10`="Intermediate tumor",
                                  `15`="Vascular tumor",
                                  `7`="Mixed tumor 1",
                                  `9`="Mixed tumor 2", 
                                  `6`="MPO+ infiltrates", 
                                  `12`="MES-1 enriched")
integrated_seurat$neighborhood_names <- Idents(integrated_seurat)

for_plot <- prop.table(table(HGG.sce$cn_celltypes1, HGG.sce$annotations_lv2), 
                       margin = 1)

pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")
pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "row")

####Hypergeometric enrichment test 
neighborhoods <- unique(integrated_seurat$neighborhoods)
annotations <- unique(integrated_seurat$annotations_lv2)
p_value_matrix <- matrix(NA, nrow = length(neighborhoods), ncol = length(annotations),
                         dimnames = list(neighborhoods, annotations))
n = ncol(integrated_seurat)
for (i in 1:length(neighborhoods)) {
  print(i)
  for (j in 1:length(annotations)) {
    print(j)
    # Get counts for each combination of categories
    obs_count <- sum(integrated_seurat$neighborhoods == neighborhoods[i] & 
                       integrated_seurat$annotations_lv2 == annotations[j])
    
    # Perform hypergeometric test
    p_value <- phyper(obs_count - 1, 
                      sum(integrated_seurat$annotations_lv2 == annotations[j]),
                      n-sum(integrated_seurat$annotations_lv2 == annotations[j]),
                      sum(integrated_seurat$neighborhoods == neighborhoods[i]),
                      lower.tail = FALSE)
    print(p_value)
    p_value_matrix[i, j] <- p_value
  }
}
write.table(p_value_matrix, file = "HG_Pvalue.txt", sep = '\t', quote = F)
long_df <- melt(p_value_matrix, variable.name = "Variable", value.name = "Value")
colnames(long_df) <- c("CN", "Cell", "pval")
long_df$fdr = p.adjust(long_df$pval, method = "BH")


######
#Generate Neighborhood Masks
######
annotations_col <- "neighborhoods"
cluster_factor <- as.numeric(as.factor(integrated_seurat$neighborhoods))
integrated_seurat$cluster_levels <- cluster_factor
conversion_table <- data.frame(
  Level = levels(as.factor(integrated_seurat$neighborhoods)),
  Numeric = levels(factor(as.numeric(factor(cluster_factor))))
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
  write.csv(annotations_df, file = paste0("Annotation_Masks/Filtered_Neighborhood",name,"_",annotations_col,"_",pix,"_11172023.csv"))
}


