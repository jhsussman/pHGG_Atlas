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
options(future.globals.maxSize = 80000 * 1024^2)
setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/Multiome")

whole_hgg_obj_093022 <- readRDS("Reference_RNA_Original.rds")
multiome_data = readRDS("Seurat_Multiome_All_DoubletsRemoved.rds")

DefaultAssay(multiome_data) = "RNA"
multiome_data[["SCT"]] = NULL
multiome_data = JoinLayers(multiome_data)
multiome_data <- SCTransform(multiome_data, verbose = TRUE)

anchors <- FindTransferAnchors(
  reference = whole_hgg_obj_093022,
  query = multiome_data,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50
)
mapped_data <- MapQuery(
  anchorset = anchors,
  query = multiome_data,
  reference = whole_hgg_obj_093022,
  refdata = list(
    celltype.proj = "cellType"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap2"
)
DimPlot(mapped_data, group.by = "predicted.celltype.proj", 
        reduction = 'ref.umap', raster = F) + coord_fixed()
table(mapped_data$predicted.celltype.proj)

#All merged
DefaultAssay(whole_hgg_obj_093022) <- "RNA"
DefaultAssay(mapped_data) <- "RNA"
whole_hgg_obj_093022$id <- "reference"
whole_hgg_obj_093022$predicted.celltype.proj <- 'reference'
mapped_data$id <- "query"
refquery <- merge(whole_hgg_obj_093022, mapped_data)
refquery[["pca"]] <- merge(whole_hgg_obj_093022[["pca"]], mapped_data[["ref.pca"]])
refquery[["umap"]] <- merge(whole_hgg_obj_093022[["umap2"]], mapped_data[["ref.umap"]])
refquery$merged_cellType <- c(whole_hgg_obj_093022$cellType, mapped_data$predicted.celltype.proj)
refquery$datatime <- c(rep("Initial", length(whole_hgg_obj_093022$cellType)), 
                       rep("Added", length(mapped_data$predicted.celltype.proj)))
Idents(refquery) = "merged_cellType"
refquery <- RenameIdents(refquery, "Endothelial Cell"="Endothelial Cells", 
                        "Macrophage/Microglia"="Macrophages/Microglia",
                        "Mature Neuron"="Mature Neurons",
                        "Mature Oligodendrocyte"="Mature Oligodendrocytes",
                        "Neuroglial"="Other Neural + Glial",
                        "Pericyte"="Pericytes", 
                        "T Cell"="T Cells")
refquery$merged_cellType <- Idents(refquery)

myColors = c("Endothelial Cells" = "#94CB5E", "Macrophages/Microglia" = "#048CFC", 
                      "Mature Neurons" = "#FC948C", "Mature Oligodendrocytes" = "#B4C4FC", 
                      "Other Neural + Glial" = "#AB130C", "Pericytes" = "#FCC404", "T Cells"= "#AC04B4")
myColors['reference'] = '#F9F9F9'
myColors_1<- myColors

refquery$predicted.celltype.proj = factor(refquery$predicted.celltype.proj)
refquery[["SCT"]] = NULL
DefaultAssay(refquery) = "RNA"
refquery[["integrated"]] = NULL
refquery$predicted.celltype.proj = factor(refquery$predicted.celltype.proj, 
                                        levels = c("reference", "GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                                                   "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like"))
p1 = DimPlot(refquery, group.by = 'predicted.celltype.proj', cols = myColors_1, raster = F, reduction = "umap") + coord_fixed() 
p1


table(refquery$timepoint)
DimPlot(refquery, reduction = 'umap', raster = F, 
        group.by = 'merged_cellType', alpha = 0.1) + coord_fixed()
p1 = FeaturePlot(mapped_data, 
            features = 'predicted.celltype.proj.score', order = T,
            raster = TRUE, raster.dpi = c(3000,3000), pt.size = 5, alpha = 1) + coord_fixed() +
  scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdBu"))) +theme_void()
ggsave(p1, filename="Figures/Multiome_Projection_Confidence.pdf", device="pdf", 
       width = 20, height = 20, units = "cm")  

p2 = DimPlot(whole_hgg_obj_093022, reduction = "umap2", 
                 group.by = 'id', order = T, cols = '#F9F9F9', 
                 raster = TRUE, raster.dpi = c(3000,3000), pt.size = 5, alpha = 1) + coord_fixed() + theme_void()
ggsave(p2, filename="Figures/Multiome_Projection_Confidence_Background.pdf", device="pdf", 
       width = 20, height = 20, units = "cm")  

#saveRDS(mapped_data, "Seurat_Multiome_Annotated.rds")

#Run copy number variations 
mapped_data = readRDS("Seurat_Multiome_Annotated.rds")
table(mapped_data$predicted.celltype.proj)
Idents(mapped_data) = "predicted.celltype.proj"
mapped_data = RenameIdents(mapped_data, "Endothelial Cell" = "Vascular", 
                           "Pericyte" = "Vascular", "Macrophage/Microglia" = "White Blood Cells", 
                           "T Cell" = "White Blood Cells", "Mature Oligodendrocyte" = "Mature Glia")
mapped_data$For_inferCNV = Idents(mapped_data)

gof <- read.table("hg38_gencode_v27.txt")
options(scipen = 100) 

rawCounts = mapped_data[["RNA"]]$counts
annotations <- mapped_data@meta.data
annotations$cellnames <- colnames(mapped_data)
annotations <- annotations %>% dplyr::select(For_inferCNV)
annotations <- as.matrix(annotations, row.names = FALSE)
inferCNV_object <- CreateInfercnvObject(raw_counts_matrix = rawCounts, 
                                    annotations_file = annotations, 
                                    delim="\t", 
                                    gene_order_file= "hg38_gencode_v27.txt", 
                                    ref_group_names= c("Mature Glia", "Vascular", "White Blood Cells"))
inferCNV_object <- infercnv::run(inferCNV_object, cutoff=0.1, 
                              out_dir= "InferCNV", 
                              cluster_by_groups = F, denoise=T, HMM=T, 
                              leiden_resolution = 0.0005, 
                              analysis_mode = "subclusters",
                              output_format = "png", 
                              num_threads = 24)


