source('scDataAnalysis_Utilities_simp.R')
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# get gene annotations for hg38
#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotation) <- "UCSC"
#genome(annotation) <- "hg38"


dir0 = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/pHGG/final_output/'
samples = dir(dir0)
samples = samples[grepl(samples, pattern = '^HTAN_pHGG')]
annotation = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/scGEL/data/annotation_hg38.rds')
seu_list = list()
for(sample0 in samples){
  mfile = paste0(dir0, sample0, '/filtered_matrix/MACS2/FILTER/reConstruct_matrix/matrix.rds')
  if(!file.exists(mfile)) next
  tmp = readRDS(mfile)
  cs = colSums(tmp > 0)
  tmp = tmp[, cs > 2000]
  rs = rowSums(tmp > 0)
  tmp = tmp[rs > 0, ]
  
  # add qc metric
  qc_stat_file = paste0(dir0, sample0, '/summary/', sample0, '.MACS2.qc_per_barcode.txt')
  if(file.exists(qc_stat_file)){
    qc_singlecell = fread(qc_stat_file)
    qc_singlecell = qc_singlecell[bc %in% colnames(tmp)]
    qc_singlecell = data.frame(qc_singlecell)
    rownames(qc_singlecell) = qc_singlecell$bc
    qc_singlecell$bc = NULL
    names(qc_singlecell) =  c("total.unique.frags", "frac.mito",  "frac.peak",
                              "frac.promoter", "frac.tss", "frac.enhancer", "tss_enrich_score")
    #seurat.obj <- AddMetaData(seurat.obj, metadata = qc_singlecell)
  }
  
  chrom_assay <- CreateChromatinAssay(
    counts = tmp,
    fragments = paste0(dir0, sample0, '/summary/', sample0, '.fragments.tsv.gz'),
    annotation = annotation
  )
  # create a seurat obj for each sample
  seurat.obj = CreateSeuratObject(counts = chrom_assay, assay = 'ATAC', 
                                  meta.data = qc_singlecell)
  
  seurat.obj$sampleName = gsub('HTAN_pHGG_', '', sample0)
  seu_list[[sample0]] = seurat.obj
}

samples = names(seu_list)
seurat.atac0 = merge(seu_list[[1]], seu_list[-1],
                     add.cell.ids = sapply(samples, function(x) gsub('HTAN_pHGG_', '', x)))

## add more metadata
sample.sheet = subset(seurat.atac0@meta.data, select = c(sampleName)) %>% unique()
sample.sheet = data.table(sample.sheet)
sample.sheet[, 'sampleID' := unlist(strsplit(sampleName, '_'))[1], by = sampleName]
sample.sheet[, 'Region' := unlist(strsplit(sampleName, '_'))[2], by = sampleName]

file.ann <- openxlsx::read.xlsx('MetaData/HTAN_pHGG_file_annot.xlsx', sheet = 1)
file.ann = data.table(file.ann)
file.ann[, 'sampleID' := gsub('7316-', '', sampleID)]
file.ann = subset(file.ann, select = c('sampleID', 'patientID', 'timepoint', 'molecularClass')) %>% unique()
setkey(file.ann, sampleID)
setkey(sample.sheet, sampleName)

sample.sheet$timepoint = file.ann[sample.sheet$sampleID]$timepoint
sample.sheet$patientID = file.ann[sample.sheet$sampleID]$patientID

seurat.atac0$timepoint = sample.sheet[seurat.atac0$sampleName]$timepoint
seurat.atac0$sampleID = sample.sheet[seurat.atac0$sampleName]$sampleID
seurat.atac0$patientID = sample.sheet[seurat.atac0$sampleName]$patientID
seurat.atac0$region = sample.sheet[seurat.atac0$sampleName]$Region


DefaultAssay(seurat.atac0) = 'ATAC'

saveRDS(seurat.atac0, file = 'SeuratObjects/seurat_atac_withFragments.rds')
sampleISeDs = seurat.atac0$sampleID
mtx = seurat.atac0@assays$ATAC@counts

## remove peaks in chrM
rnames = rownames(mtx)
pks.nochrm = rnames[!grepl(rnames, pattern = '^chrM')]
mtx = mtx[pks.nochrm, ]


## alternatively
seurat.atac0 = readRDS('SeuratObjects/seurat_atac_withFragments.rds')
mtx.ga <- GeneActivity(seurat.atac0)
saveRDS(mtx.ga, file = 'SeuratObjects/mtx_gene_activity_signac.rds')


freq_by_sample <- sapply(unique(sampleIDs), function(x){
  return(Matrix::rowMeans(mtx[, sampleIDs == x] > 0))
})
freq.max = apply(freq_by_sample, 1, max)
names(freq.max) = rownames(freq_by_sample)
filtered.pks1 = names(which(freq.max < 0.01))
mtx = mtx[!rownames(mtx) %in% filtered.pks1, ]

npc = 50
rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.002))

nvap = 10000
inherited.mdata <- seurat.atac0@meta.data
seurat.atac = doBasicSeurat_atac_updated(mtx, npc = npc, 
                                         norm_by = 'tf-idf',
                                         top.variable = nvap,
                                         regressOnPca = TRUE,
                                         reg.var = 'nFeature_ATAC', 
                                         excludePks.fromVAP = filtered.pks, 
                                         meta.data =inherited.mdata)

vaps1 = VariableFeatures(seurat.atac)

seurat.atac <- RunUMAP(seurat.atac, dims = 1:npc)

DimPlot(seurat.atac, group.by = 'sampleID', label = T) 
DimPlot(seurat.atac, group.by = 'timepoint', label = T) 

saveRDS(seurat.atac, file = 'SeuratObjects/seurat_atac_pooled.rds')

## reselect variable features across clusters ####
npc = 50
seurat.atac = readRDS('SeuratObjects/seurat_atac_pooled.rds')
seurat.atac <- FindNeighbors(seurat.atac, reduction = 'pca', 
                             dims = 1:npc)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.2)

nvap = 10000
niter = 1
k = 0 
repeat{
  k = k + 1
  if(k > niter) break
  mtx_by_cls <- sapply(unique(seurat.atac$seurat_clusters), function(x) {
    cl_data <- 1* (mtx[, seurat.atac$seurat_clusters == x] > 0)
    return(Matrix::rowSums(cl_data))
  })
  
  mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
  
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sds = sort(sds, decreasing = T)
  sele.features = names(sds[1:nvap])
  sele.features = setdiff(sele.features, filtered.pks)
  #mtx.norm = TF_IDF(mtx[sele.features, ])
  
  #tmp <- mtx[setdiff(rownames(mtx), sele.features), ]
  #data0 <- rbind(mtx.norm, tmp)
  #seurat.atac[['ATAC']]@data = data0[rownames(mtx), ]
  
  VariableFeatures(seurat.atac) = sele.features
  seurat.atac <- ScaleData(seurat.atac, features = sele.features)
  seurat.atac <- RunPCA(seurat.atac, npcs = npc, 
                        features = VariableFeatures(seurat.atac))
  seurat.atac <- regress_on_pca(seurat.atac, reg.var = 'nFeature_ATAC')
  seurat.atac <- RunUMAP(seurat.atac, dims = 1:npc)
  
  
  seurat.atac <- FindNeighbors(seurat.atac, reduction = 'pca', dims = 1:npc)
  seurat.atac <- FindClusters(seurat.atac, resolution = 0.2)
  
}

DimPlot(seurat.atac, group.by = 'sampleID', label = T) 
DimPlot(seurat.atac, group.by = 'timepoint', label = T) 



## do harmony integration
library(harmony)
seurat.atac <- RunHarmony(seurat.atac, 'sampleID', assay.use = 'ATAC')
seurat.atac <- RunUMAP(seurat.atac, reduction = 'harmony', 
                       reduction.name = 'harmonyUMAP',  reduction.key = 'harmonyUMAP_', dims = 1:30)
DimPlot(seurat.atac, group.by = 'sampleID', label = T, repel = T, reduction = 'harmonyUMAP')
DimPlot(seurat.atac, group.by = 'timepoint', label = T, repel = T, reduction = 'harmonyUMAP')

seurat.atac <- FindNeighbors(seurat.atac, reduction = 'harmony', dims = 1:30, 
                             graph.name = c('harmony_nn', 'harmony_snn'))
seurat.atac <- FindClusters(seurat.atac, resolution = 0.1, graph.name = 'harmony_snn')
seurat.atac <- FindClusters(seurat.atac, resolution = 0.2, graph.name = 'harmony_snn')

niter = 1
k = 0 
repeat{
  k = k + 1
  if(k > niter) break
  mtx_by_cls <- sapply(unique(seurat.atac$seurat_clusters), function(x) {
    cl_data <- 1* (mtx[, seurat.atac$seurat_clusters == x] > 0)
    return(Matrix::rowSums(cl_data))
  })
  
  mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
  
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sds = sort(sds, decreasing = T)
  sele.features = names(sds[1:nvap])
  sele.features = setdiff(sele.features, filtered.pks)
  #mtx.norm = TF_IDF(mtx[sele.features, ])
  
  #tmp <- mtx[setdiff(rownames(mtx), sele.features), ]
  #data0 <- rbind(mtx.norm, tmp)
  #seurat.atac[['ATAC']]@data = data0[rownames(mtx), ]
  
  VariableFeatures(seurat.atac) = sele.features
  seurat.atac <- ScaleData(seurat.atac, features = sele.features)
  seurat.atac <- RunPCA(seurat.atac, npcs = npc, 
                        features = VariableFeatures(seurat.atac))
  seurat.atac <- regress_on_pca(seurat.atac, reg.var = 'nFeature_ATAC')
  seurat.atac <- RunUMAP(seurat.atac, dims = 1:npc)
  
  
  seurat.atac <- FindNeighbors(seurat.atac, reduction = 'pca', dims = 1:npc)
  seurat.atac <- FindClusters(seurat.atac, resolution = 0.2)
  
}

seurat.atac <- RunHarmony(seurat.atac, 'sampleID', assay.use = 'ATAC')
seurat.atac <- RunUMAP(seurat.atac, reduction = 'harmony', 
                       reduction.name = 'harmonyUMAP',  reduction.key = 'harmonyUMAP_', dims = 1:30)
DimPlot(seurat.atac, group.by = 'sampleID', label = T,  reduction = 'harmonyUMAP')
DimPlot(seurat.atac, group.by = 'timepoint', label = T,  reduction = 'harmonyUMAP')

seurat.atac <- FindNeighbors(seurat.atac, reduction = 'harmony', dims = 1:30, 
                             graph.name = c('harmony_nn', 'harmony_snn'))
seurat.atac <- FindClusters(seurat.atac, resolution = 0.1, graph.name = 'harmony_snn')
seurat.atac <- FindClusters(seurat.atac, resolution = 0.2, graph.name = 'harmony_snn')

## add gene activity score
seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = mtx.ga)

DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))
seurat.atac <- FindVariableFeatures(seurat.atac)
DefaultAssay(seurat.atac) <- 'ATAC'

saveRDS(seurat.atac, paste0('SeuratObjects/seurat_atac_vap', nvap, '_filteredLowQCells.rds'))


