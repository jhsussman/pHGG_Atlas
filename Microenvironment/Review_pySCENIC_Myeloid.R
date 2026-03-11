library(Seurat)
library(tidyverse)
library(data.table)
library(SCopeLoomR)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(ggridges)
library(heatmaply)
library(AUCell)
library(SCENIC)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(Matrix)

setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Macrophages/SCENIC")

out_dir <- "pySCENIC_Output"
out_files <- list.files(out_dir, full.names = TRUE)
seurat_path <- "../Reanalysis_scRNA_macrophage_filtered2.rds"
analysis <- "pHGG_Myeloid"
path2loom <- out_files[grep("pySCENIC_AUCell.loom", out_files)]
loom <- open_loom(path2loom, mode = "r")
regulonsMat <- get_regulons(loom, column.attr.name = "Regulons")
dim(regulonsMat)
write.csv(regulonsMat,file = paste0(out_dir, "/", analysis, "_pySCENIC_regulons_incidence_matrix.csv"))
regulons <- SCENIC::regulonsToGeneLists(regulonsMat) 
head(regulons)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
saveRDS(regulonsAUC, file = "regulonsAUC.rds")
close_loom(loom)

### Add desired data to original Seurat object (RNA only)
myeloid_rna <- readRDS(seurat_path)
AUCmat <- AUCell::getAUC(regulonsAUC)
AUCmat[1:10, 1:10]
rownames(AUCmat)
myeloid_rna[['SCENIC']] <- CreateAssayObject(data = AUCmat)
DefaultAssay(myeloid_rna) <- "SCENIC"
myeloid_rna <- ScaleData(myeloid_rna, assay = 'SCENIC', features = rownames(AUCmat))
#saveRDS(myeloid_rna, file = paste0(out_dir, "/", analysis, "_SCENIC.rds"))

###
myeloid_rna <- readRDS("pySCENIC_Output/pHGG_Myeloid_SCENIC.rds")

Idents(myeloid_rna) <- "cell_labels1"
deg <- FindAllMarkers(myeloid_rna, logfc.threshold = 0.01, only.pos = TRUE)
write.table(deg, file = "SCENIC_Differential_Scores.tsv", sep = '\t', quote = F)
deg %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) -> top
top = top$gene
DoHeatmap(myeloid_rna, features = as.character(unique(top$gene)), group.by = 'cell_labels1')+
  scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave(p1, filename = "SCENIC_Heatmap.png", device="png", width = 20, height = 12)

myeloid_rna@assays$SCENIC@counts = myeloid_rna@assays$SCENIC@data
combined_averages <- AverageExpression(myeloid_rna, return.seurat = FALSE, assays = "SCENIC", group.by = "cell_labels1")
combined_averages <- data.frame(as(combined_averages$SCENIC, "matrix"))
subset_matrix <- combined_averages[top, ]
pheatmap(
  subset_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = FALSE, 
  cluster_cols = FALSE, fontsize_row = 8.5,
  labels_row = top, scale = "row",
  main = "Average Regulon AUC"  
)