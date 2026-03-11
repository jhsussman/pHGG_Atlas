library(ggplot2)
library(Seurat)
library(Signac)
library(ggrastr)
library(dittoSeq)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(paletteer)
library(monocle3)
library(SeuratWrappers)

tumor_colors <- paletteer_d("ggthemes::Classic_10", n = 10)
tumor_colors <- as.vector(tumor_colors)
names(tumor_colors) <- c("GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                         "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")

setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Analysis_Revisions")

ng_rna <- readRDS("../snRNA-seq/seurat_rna_Neuroglia_final.rds")
ng_rna <- AddMetaData(ng_rna, metadata=readRDS("../snRNA-seq/metadata_seurat_rna_Neuroglia_final.rds"))
Idents(ng_rna) <- "timepoint"
ng_rna <- RenameIdents(ng_rna, "Initial CNS Tumor" = "Initial Tumor", "Recurrence" = "Timepoint_2", "Progressive (Non-Autopsy)" = "Timepoint_2", 
                       "Progressive (Autopsy)" = "Timepoint_2")
ng_rna$timemerge <- Idents(ng_rna)
Idents(ng_rna) <- "timepoint"
ng_rna <- RenameIdents(ng_rna, "Initial CNS Tumor" = "Initial Tumor", "Recurrence" = "Recurrence/Progression", "Progressive (Non-Autopsy)" = "Recurrence/Progression", 
                       "Progressive (Autopsy)" = "Autopsy")
ng_rna$threetime <- Idents(ng_rna)

#Monocle3 
cds <- as.cell_data_set(ng_rna)
cds <- cluster_cells(cds)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,cds@colData$cell_state0 == "AC-like 2"]))
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cell_state1",
           label_cell_groups = TRUE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60") + coord_fixed()
#saveRDS(cds, "CDS_Neuroglia_Monocle3.rds")

cds = readRDS("CDS_Neuroglia_Monocle3.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = TRUE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60") + coord_fixed()

new_metadata=readRDS("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Analysis_Revisions/Metadata_cell_state1.rds")
integrated.sub <- as.Seurat(cds, assay = NULL)
integrated.sub = AddMetaData(integrated.sub, metadata = new_metadata)
FeaturePlot(integrated.sub, "monocle3_pseudotime") + coord_fixed()

integrated.sub$cell_state1 = factor(integrated.sub$cell_state1, levels = rev(c("GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                                                                               "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")))
p1 = RidgePlot(integrated.sub, features = "monocle3_pseudotime",group.by = "cell_state1", cols = tumor_colors) + theme_bw()
ggsave(p1, filename = "Figures/Pseudotime_Ridgeplot.pdf", width = 6.5, height = 3)
