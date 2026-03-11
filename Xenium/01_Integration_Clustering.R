source("/mnt/isilon/tan_lab/sussmanj/Xenium/Xenium_Tools.R")

setwd("/mnt/isilon/tan_lab/sussmanj/Xenium/pHGG")
options(future.globals.maxSize = 80000 * 1024^2)
options(Seurat.object.assay.version = "v5")

#Fix the xenium output
prefix = "/mnt/isilon/tan_lab/sussmanj/Xenium/pHGG/"
xenium.dir = paste0(prefix, "output-XETG00217__0033486__REG2__20240912__195600")
transcripts <- arrow::read_parquet(file.path(xenium.dir, "transcripts.parquet"))
gz_con <- gzfile(file.path(xenium.dir, "transcripts.csv.gz"), "wt")
write.csv(transcripts, gz_con, row.names = FALSE, quote = TRUE)
close(gz_con)
xenium.dir = paste0(prefix, "output-XETG00217__0033486__REG1__20240912__195600")
transcripts <- arrow::read_parquet(file.path(xenium.dir, "transcripts.parquet"))
gz_con <- gzfile(file.path(xenium.dir, "transcripts.csv.gz"), "wt")
write.csv(transcripts, gz_con, row.names = FALSE, quote = TRUE)
close(gz_con)
xenium.dir = paste0(prefix, "output-XETG00217__0033487__Region_1__20240912__195600")
transcripts <- arrow::read_parquet(file.path(xenium.dir, "transcripts.parquet"))
gz_con <- gzfile(file.path(xenium.dir, "transcripts.csv.gz"), "wt")
write.csv(transcripts, gz_con, row.names = FALSE, quote = TRUE)
close(gz_con)
xenium.dir = paste0(prefix, "output-XETG00217__0033487__Region_2__20240912__195600")
transcripts <- arrow::read_parquet(file.path(xenium.dir, "transcripts.parquet"))
gz_con <- gzfile(file.path(xenium.dir, "transcripts.csv.gz"), "wt")
write.csv(transcripts, gz_con, row.names = FALSE, quote = TRUE)
close(gz_con) 

#Load xenium
prefix = "/mnt/isilon/tan_lab/sussmanj/Xenium/pHGG/"
xenium.dir = paste0(prefix, "output-XETG00217__0033487__Region_1__20240912__195600")
pHGG_161 = LoadXenium(data.dir = xenium.dir) 
pHGG_161$orig.ident = "pHGG_161"
xenium.dir = paste0(prefix, "output-XETG00217__0033487__Region_2__20240912__195600")
pHGG_6477 = LoadXenium(data.dir = xenium.dir) 
pHGG_6477$orig.ident = "pHGG_6477"
xenium.dir = paste0(prefix, "output-XETG00217__0033486__REG2__20240912__195600")
pHGG_339 = LoadXenium(data.dir = xenium.dir) 
pHGG_339$orig.ident = "pHGG_339"
xenium.dir = paste0(prefix, "output-XETG00217__0033486__REG1__20240912__195600")
pHGG_942 = LoadXenium(data.dir = xenium.dir) 
pHGG_942$orig.ident = "pHGG_942"

seurat_list = list(pHGG_161, pHGG_6477, pHGG_339, pHGG_942)
#saveRDS(seurat_list, "Seurat_Objects/Seurat_List_4.rds")

seurat_list = readRDS("Seurat_Objects/Seurat_List_4.rds")
pHGG_161 = seurat_list[[1]]
pHGG_6477 = seurat_list[[2]]
pHGG_339 = seurat_list[[3]]
pHGG_942 = seurat_list[[4]]

pHGG_161[["BlankCodeword"]] = NULL; pHGG_161[["ControlCodeword"]] = NULL; pHGG_161[["ControlProbe"]] = NULL
pHGG_6477[["BlankCodeword"]] = NULL; pHGG_6477[["ControlCodeword"]] = NULL; pHGG_6477[["ControlProbe"]] = NULL
pHGG_339[["BlankCodeword"]] = NULL; pHGG_339[["ControlCodeword"]] = NULL; pHGG_339[["ControlProbe"]] = NULL
pHGG_942[["BlankCodeword"]] = NULL; pHGG_942[["ControlCodeword"]] = NULL; pHGG_942[["ControlProbe"]] = NULL

plot_points(pHGG_161, alpha = 0.05)
plot_points(pHGG_6477, alpha = 0.1)
plot_points(pHGG_339, alpha = 0.1)
plot_points(pHGG_942, alpha = 0.1)

#Crop out images 
pHGG_161 = filter_xenium_seurat(xenium_seurat = pHGG_161, coords_csv = "pHGG_161_ROI.csv")
pHGG_942 = filter_xenium_seurat(xenium_seurat = pHGG_942, coords_csv = "pHGG_942_ROI.csv")
pHGG_339 = filter_xenium_seurat(xenium_seurat = pHGG_339, coords_csv = "pHGG_339_ROI.csv")
pHGG_6477 = filter_xenium_seurat(xenium_seurat = pHGG_6477, coords_csv = "pHGG_6477_ROI.csv")

plot_points(pHGG_161, alpha = 0.05)
plot_points(pHGG_6477, alpha = 0.1)
plot_points(pHGG_339, alpha = 0.1)
plot_points(pHGG_942, alpha = 0.5)

seurat_list = list(pHGG_161, pHGG_6477, pHGG_339, pHGG_942)
seurat_merged = merge(seurat_list[[1]], y = seurat_list[-1])
saveRDS(seurat_merged, "Seurat_Objects2/Seurat_Merged_4_Cropped.rds")

VlnPlot(seurat_merged, group.by = "orig.ident", features = c("nCount_Xenium", "nFeature_Xenium"), ncol = 2, pt.size = 0)

#Do Basic Seurat 
seurat_merged
xenium.obj <- subset(seurat_merged, subset = nCount_Xenium > 10 & nFeature_Xenium > 5)
xenium.obj
table(seurat_merged$orig.ident)
table(xenium.obj$orig.ident)

#Integrate using Log Normalize
xenium.obj = JoinLayers(xenium.obj)
DefaultAssay(xenium.obj) = "Xenium"
xenium.obj[["SCT"]] = NULL
xenium.obj[["Xenium"]] <- split(xenium.obj[["Xenium"]], f = xenium.obj$orig.ident)
xenium.obj <- NormalizeData(xenium.obj, assay = "Xenium")
VariableFeatures(xenium.obj) <- rownames(xenium.obj)
xenium.obj <- ScaleData(xenium.obj, vars.to.regress = "nCount_Xenium")
xenium.obj <- RunPCA(xenium.obj, features = rownames(xenium.obj))
integrated_seurat <- IntegrateLayers(object = xenium.obj, method = RPCAIntegration, features = rownames(xenium.obj), 
                                     verbose = T, orig = "pca", dims = 1:30,
                                     new.reduction = "integrated.rpca")
integrated_seurat <- JoinLayers(integrated_seurat)
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "integrated.rpca", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = c(0.025)) 
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30, reduction = "integrated.rpca", 
                             reduction.name = "umap.rpca", return.model = TRUE)
DimPlot(integrated_seurat, group.by = "orig.ident", reduction = "umap.rpca", raster = F) + coord_fixed()
DimPlot(integrated_seurat, group.by = "Xenium_snn_res.0.2", reduction = "umap.rpca", raster = F, label = T) + coord_fixed()
#saveRDS(integrated_seurat, file = "Seurat_Objects2/Merged_Seurat_Process_RPCA_C10_F5_Log.RDS")

#Plotting 
cluster_counts <- table(integrated_seurat$Xenium_snn_res.0.3)
integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
  mutate(Xenium_snn_res.0.3 = ifelse(Xenium_snn_res.0.3 %in% names(cluster_counts[cluster_counts < 10]), 
                                     "Other", Xenium_snn_res.0.3))
DimPlot(integrated_seurat, group.by = "Xenium_snn_res.0.3", reduction = "umap.rpca", raster = F, label = T, label.box = T, alpha = 0.4) + coord_fixed()

FeaturePlot(integrated_seurat, features = c("nCount_Xenium", "nFeature_Xenium"), 
            reduction = "umap.rpca", raster = F, max.cutoff = 'q99', order = T) + coord_fixed()

#Explore 
integrated_seurat = readRDS("Seurat_Objects/Merged_Seurat_Process_RPCA_C10_F5_Log.RDS")
pdf("Xenium_Marker_Plots_Log.pdf")
for(marker in rownames(integrated_seurat)){
  print(marker)
  try(p <- FeaturePlot(integrated_seurat, features = marker, raster = T, reduction = "umap.rpca", raster.dpi = c(1200, 1200),
                       pt.size = 1.5, alpha = 1, order = T)) + coord_fixed()
  try(print(p))
}
dev.off()

rownames(integrated_seurat)
FeaturePlot(integrated_seurat, features = c("nCount_Xenium", "nFeature_Xenium"), reduction = "umap.rpca", max.cutoff = 'q99', ncol = 1) 

#Differential expression 
integrated_seurat = readRDS("Merged_Seurat_Process_RPCA_nCount10.RDS")

DefaultAssay(integrated_seurat) <- "Xenium"
integrated_seurat <- NormalizeData(integrated_seurat)
integrated_seurat <- JoinLayers(integrated_seurat)
Idents(integrated_seurat) = "SCT_snn_res.0.3"
markers = FindAllMarkers(integrated_seurat, only.pos = T)
markers$diff.pct = markers$pct.1 - markers$pct.2

top <- markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
DotPlot(integrated_seurat, features = unique(top$gene), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


