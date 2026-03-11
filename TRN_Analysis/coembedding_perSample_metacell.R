## transfer label from scRNA to scATAC using Seurat 

source('scDataAnalysis_Utilities_simp.R')
library(gridExtra)
`%notin%` = Negate(`%in%`)

args = commandArgs(T)
sampleID0 = args[1]

seuratAtacPath = 'SeuratObjects/seurat_atac_signacIntegrated_bySample_withReference_4036_4037.rds'
reduction2use = 'pca' # pca or harmony

## load scRNA data and cell type annotation ####
seurat.rna <- readRDS('SeuratObjects/seurat_rna_chop.rds')

seurat.rna = subset(seurat.rna, sampleID == paste0('7316-', sampleID0))

seurat.atac = readRDS(seuratAtacPath)
seurat.atac = subset(seurat.atac, sampleID == sampleID0)

## recreate seurat object using the sample only
mtx = seurat.atac@assays$ATAC@counts
npc = 30
rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.003))

nvap = 10000
inherited.mdata <- seurat.atac@meta.data
seurat.atac = doBasicSeurat_atac_updated(mtx, npc = npc, 
                                         norm_by = 'logNormalize',
                                         top.variable = nvap,
                                         regressOnPca = TRUE,
                                         reg.var = 'nFeature_ATAC', 
                                         excludePks.fromVAP = filtered.pks, 
                                         meta.data =inherited.mdata)



## use GAS = promote + gene body accessibility for label transfer ####
if(all(names(seurat.atac@assays) != 'ACTIVITY')){
  mtx.ga <- readRDS('SeuratObjects/mtx_gene_activity_signac.rds')
  mtx.ga <- mtx.ga[, colnames(seurat.atac)]
  
  seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = mtx.ga)
  seurat.atac <- NormalizeData(
    object = seurat.atac,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.atac$nCount_ACTIVITY)
  )
  
}
DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

## transfer label 
DefaultAssay(seurat.rna) = 'RNA'
#seurat.rna = LogNormalize(seurat.rna)
seurat.rna = FindVariableFeatures(seurat.rna, nfeatures = 2200)
vgenes = VariableFeatures(seurat.rna)
gfreq = rowMeans(seurat.rna@assays$RNA@counts > 0)
rgenes = names(which(gfreq < 0.0025))
vgenes = setdiff(vgenes, rgenes)
VariableFeatures(seurat.rna) <- vgenes
seurat.rna = ScaleData(seurat.rna) %>% RunPCA(verbose = F, npcs = 30) %>% RunUMAP(reduction = 'pca', dim = 1:30)
p0 <- DimPlot(seurat.rna, group.by = 'cellType')

genes4anchors = VariableFeatures(object = seurat.rna)

transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = "RNA",
                                        normalization.method = 'LogNormalize',
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)


celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.rna$cellType,
                                     weight.reduction = seurat.atac[['pca']],
                                     dims = 1:ncol(seurat.atac[['pca']]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, select = c('predicted.id', 'prediction.score.max'))
names(celltype.predictions) = c('predicted.cellType.perSample', 'prediction.score.max.perSample')
seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)

DefaultAssay(seurat.atac) <- 'ATAC'
seurat.atac <- RunUMAP(seurat.atac, reduction = 'pca', dims = 1:npc)
p1 <- DimPlot(seurat.atac, group.by = 'predicted.cellType.perSample')

message('Label trasfer done! Now saving predicted data...')
saveRDS(celltype.predictions, file = paste0('MetaData/cellType_Prediction/cellType_Prediction_', sampleID0, '.rds' ))
message('Now working on co-embedding ... ')

## co-embedding ####
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <- GetAssayData(seurat.rna, assay = "RNA", 
                        slot = "data")[genes4anchors, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.atac[["pca"]],
                           dims = 1:ncol(seurat.atac[["pca"]]))

# this line adds the imputed data matrix to the seurat.atac object
seurat.atac[["RNA"]] <- imputation
coembed <- merge(x = seurat.rna, y = seurat.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE, npcs = 30)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- FindNeighbors(coembed, dims = 1:30, reduction = 'pca')
coembed <- FindClusters(coembed, resl = 0.5)

p2 <- DimPlot(coembed, group.by = 'tech', label = T) 
p2

coembed$cellType[coembed$tech == 'ATAC'] = coembed$predicted.cellType.perSample[coembed$tech == 'ATAC']

DimPlot(coembed, group.by = 'cellType', label = T) 


## 1 to 1 cell matching ####
umap_coproj = coembed@reductions$umap@cell.embeddings
ac_cells <- colnames(coembed)[coembed$tech == "ATAC"]
rna_cells <- colnames(coembed)[coembed$tech == "RNA"]
umap.rna = umap_coproj[rna_cells, ]
umap.atac = umap_coproj[ac_cells, ]
final_matching <- data.table(atac_cell = ac_cells)
final_matching$atac_cell <- as.character(final_matching$atac_cell)

dist0 <- pracma::distmat(umap.atac, umap.rna)
final_matching$rna_cell <- sapply(1:nrow(umap.atac), function(x) 
  names(which.min(dist0[x, ])))

saveRDS(final_matching, paste0("Coembed_Results/atac_rna_coembedding_cell_matching_", sampleID0,
                               ".rds"))

message('Done co-embedding, now working on smoothing ... ')

## smoothing ####
if(F){
  K = 10
  seurat.rna <- FindNeighbors(seurat.rna, k.param = K, reduction = 'pca')
  knn.mat.rna = (seurat.rna@graphs$RNA_nn > 0) * 1
  
  seurat.atac <- FindNeighbors(seurat.atac, k.param = K, reduction = 'pca')
  knn.mat.atac = (seurat.atac@graphs$ATAC_nn > 0) * 1
  
  all(rowSums(knn.mat.rna) == K)
  all(rowSums(knn.mat.atac) == K)
  
  smooth.rna = seurat.rna@assays$RNA@data %*% t(knn.mat.rna)
  smooth.atac = (seurat.atac@assays$ATAC@counts > 0 ) %*% t(knn.mat.atac)
  
  smooth.rna = smooth.rna[, final_matching$rna_cell]
  smooth.atac = smooth.atac[, final_matching$atac_cell]
  
  saveRDS(smooth.rna, file = paste0("Coembed_Results/rna_smoothed_expr_", sampleID0, ".rds"))
  saveRDS(smooth.atac, file = paste0("Coembed_Results/atac_smoothed_access_", sampleID0, ".rds"))
  
  
}


## call metacell on coembedded obj ####
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

coembed <- RenameCells(coembed,  new.names = paste0(coembed$tech, '_', colnames(coembed)))

coembed <- SetupForWGCNA(
  coembed,
  gene_select = "fraction", # the gene selection approach
  #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "test", # the name of the hdWGCNA experiment
  features = genes4anchors
)

# construct metacells  in each group
coembed <- MetacellsByGroups(
  seurat_obj = coembed,
  group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
  k = 20, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  min_cells = 50,
  reduction = 'pca',
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
coembed <- NormalizeMetacells(coembed)

metacell_obj <- GetMetacellObject(coembed)

coembed <- ScaleMetacells(coembed, features=VariableFeatures(seurat.rna))
coembed <- RunPCAMetacells(coembed, features=VariableFeatures(seurat.rna), npcs = 10)
#coembed <- RunHarmonyMetacells(coembed, group.by.vars='Sample')
coembed <- RunUMAPMetacells(coembed, reduction='pca', dims=1:10)


p3 <- DimPlotMetacells(coembed, group.by = 'seurat_clusters') + ggtitle("cluster_metacell")
p3


ggsave(plot = marrangeGrob(list(p0, p1, p2, p3), nrow=2, ncol=2), 
       width = 15, height = 12, 
       filename = paste0('Coembed_Results/coembed_', sampleID0, '.png'), 
       device = 'png')

message('All Done!')

## construct and save metacell rna and atac matrices 
merged_cells = metacell_obj$cells_merged
nrna_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'RNA_'))
natac_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'ATAC_'))
sele.metacells = merged_cells[nrna_cells >= 4 & natac_cells >= 4]
nrna_cells = nrna_cells[names(sele.metacells)]
natac_cells = natac_cells[names(sele.metacells)]

#clean metacell
if(F){
  coembed@misc$test$wgcna_metacell_obj<- coembed@misc$test$wgcna_metacell_obj[, names(sele.metacells)]
  coembed <- ScaleMetacells(coembed, features=VariableFeatures(seurat.rna))
  coembed <- RunPCAMetacells(coembed, features=VariableFeatures(seurat.rna), npcs = 10)
  #coembed <- RunHarmonyMetacells(coembed, group.by.vars='Sample')
  coembed <- RunUMAPMetacells(coembed, reduction='pca', dims=1:10)
  
  
  p3 <- DimPlotMetacells(coembed, group.by = 'seurat_clusters') + ggtitle("cluster_metacell")
  p3
}

# groups of cells to combine
mask <- sapply(seq_len(length(nrna_cells)), function(x) colnames(coembed) %in%
                 unlist(strsplit(sele.metacells[x], ',', fixed = T)))
rna.mask <- sapply(seq_len(length(nrna_cells)), function(x) mask[, x]/nrna_cells[x])
rna.mtx <- coembed@assays$RNA@data %*% rna.mask # remember atac cells in RNA assay have 0 account

atac.mask <- sapply(seq_len(length(natac_cells)), function(x) mask[, x]/natac_cells[x])
atac.mtx <- coembed@assays$ATAC@data %*% atac.mask

colnames(rna.mtx) = colnames(atac.mtx) = paste0(sampleID0, '_', names(sele.metacells))

saveRDS(rna.mtx, file = paste0('Coembed_Results/metacells/rna_metacell_mtx_', sampleID0, '.rds'))
saveRDS(atac.mtx, file = paste0('Coembed_Results/metacells/atac_metacell_mtx_', sampleID0, '.rds'))


