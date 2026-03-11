library(Seurat)
library(Signac)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(patchwork)
library(pheatmap)
library(viridis)
library(pheatmap)

`%notin%` = Negate(`%in%`)


## remove features active in less than min_frac_per_cluser in first class
## do downsample to max_cellPer_cluster
## support wilcox and t test only
## test-model: one-rest; control
runDiffMotifEnrich <- function(mtx_score, clusters, test = 'wilcox',
                               fdr = 0.01, topn = 20,
                               min_frac_per_cluster = 0.2,
                               max_cell_per_clust = 300,
                               test.mode = 'one-rest', control = NULL,
                               order_by = 'pvalue', min_mean1 = 0.01, trim_perc = 0.025){
  set.seed(2023)
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  mtx_score = mtx_score[, clusters$barcode]
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    if(test.mode == 'one-rest') {
      mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    }
    if(test.mode == 'control'){
      if(is.null(control)) stop('Should specific control cluster!')
      bc2 = clusters[cluster %in% control]$barcode
      mtx2 = mtx_score[, colnames(mtx_score) %in% bc2]
      if(cluster0 %in% control) next
    }
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ], trim = trim_perc))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ], trim = trim_perc))
    
    pvs = rep(0.5, length(features))
    
    for(x in 1:length(features)){
      a1 = mtx1[x, ]
      a2 = mtx2[x, ]
      if(mu1[x] <= 0) next
      #if(mu1[x] <= mu2[x]) next
      if(length(which(!is.na(a1))) < 2 || length(which(!is.na(a2))) < 2) next
      if(mean(a1 > 0) < min_frac_per_cluster) next
      if(!is.null(max_cell_per_clust)){
        if(length(a1) > max_cell_per_clust) a1 = sample(a1, max_cell_per_clust)
        if(length(a2) > max_cell_per_clust) a2 = sample(a2, max_cell_per_clust)
      }
      
      if(test == 'wilcox') pvs[x] = wilcox.test(a1, a2, alternative = 'greater')$p.value
      if(test == 't') pvs[x] = t.test(a1, a2, alternative = 'greater')$p.value
      
    }
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster1' = cluster0,
                      'mean1' = mu1, 'mean0' = mu2, 
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    res0 = res0[mean1 > min_mean1]
    if(order_by == 'pvalue') res0 = res0[order(pv_adjust), ]
    if(order_by == 'mean') res0 = res0[order(-mean1), ]
    if(order_by == 'diff') res0 = res0[order(mean0 - mean1), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}



pearson_residual <- function(mat){
  N = sum(mat)
  rs = rowSums(mat)
  cs = colSums(mat)
  mu = rs %o% cs/N
  sadj = (1 - rs/N) %o% (1 - cs/N) ## standardize residual
  z = (mat - mu)/sqrt(mu * sadj)
  rnull = names(which(rs == 0))
  z[rnull, ] <- 0
  #z = Matrix(z, sparse = T)
  return(z)
}


## examin new label transfer restult ####
seurat.atac = readRDS('SeuratObjects/seurat_atac_tumor_signac_RemovedCls2_4.rds')
meta_label = readRDS('MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = meta_label)
p0 <- FeaturePlot(seurat.atac, features = 'prediction.state.score.max1', max.cutoff = 0.9)+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p1 <- DimPlot(seurat.atac, group.by = 'predicted.state1', label = T)

ggsave(p0 | p1, file = 'Figures/ATAC/States/Revision/label_transfer.pdf',
       width = 14, height = 6, device = 'pdf')
DefaultAssay(seurat.atac) = 'ATAC'

DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.75', label = T) # old clusters

seurat.atac = FindNeighbors(seurat.atac, reduction = 'integrated_lsi', dims = 2:30)
seurat.atac = FindClusters(seurat.atac, resolution = 0.6)
seurat.atac = FindClusters(seurat.atac, resolution = 0.8)
p2 <- DimPlot(seurat.atac, label = T)

## check enrichment of cells with high confident score
Idents(seurat.atac) = seurat.atac$ATAC_snn_res.0.8

seurat.atac0 = subset(seurat.atac, prediction.state.score.max > 0.6)

cluster_label <- table(seurat.atac0$predicted.state1, seurat.atac0$seurat_clusters)
cluster_label_norm <- pearson_residual(cluster_label)
pvs = pnorm(cluster_label_norm, lower.tail = F)
plot_mat <- -log10(pvs)
plot_mat[is.na(plot_mat)] = 0
plot_mat[is.infinite(plot_mat)] = max(plot_mat[!is.infinite(plot_mat)]) + 10
plot_mat_norm <- plot_mat %*% diag(1/colSums(plot_mat + 0.01)) 
colnames(plot_mat_norm) = colnames(plot_mat)

pc0 <- pheatmap(plot_mat_norm)

cluster_label_norm[cluster_label_norm > 30] = 30
pc1 <- pheatmap(cluster_label_norm)

ggsave(pc0, file = 'Figures/ATAC/States/Revision/enrichment_highScoreCells_normalized.pdf',
       width = 7, height = 6, device = 'pdf')
ggsave(pc1, file = 'Figures/ATAC/States/Revision/enrichment_highScoreCells.pdf',
       width = 7, height = 6, device = 'pdf')

pdf(file = 'Figures/ATAC/States/Revision/frac_highScoreCells.pdf')
aa = table(seurat.atac0$seurat_clusters)/table(seurat.atac$seurat_clusters)
barplot(sort(aa), col = 'lightblue')
abline(h = 0.05, lty = 2)
dev.off()

## assign based on enrichment of high confident labeled cells
seurat.atac$cell_state1 = ''
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('4', '14')] = 'AC-like'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c( '0')] = 'GPC-like'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('3', '7')] = 'Transition 1'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('10', '13')] = 'MES-like'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('9', '11')] = 'NEU-like'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('2')] = 'OPC/NPC-like'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('6', '12')] = 'OC-like'  
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('1')] = 'Transition 2'
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('5')] = 'HSP+ Cells'
seurat.atac$cell_state1[seurat.atac$ATAC_snn_res.0.8 %in% c('8')] = 'Unknown' 

p3 <- DimPlot(seurat.atac, group.by = 'cell_state1', label = T, label.size = 4) +
  scale_color_brewer(palette = 'Paired')
ggsave(p2 + p3, file = 'Figures/ATAC/States/Revision/cell_state1.pdf',
       width = 14, height = 6, device = 'pdf')



## update signature score using new rna cell state annotation -- use wilcox test ####
seurat.rna = readRDS('SeuratObjects/seurat_rna_Neuroglia_final.rds')
mdata = readRDS('MetaData/Metadata_cell_state1.rds')
seurat.rna = AddMetaData(seurat.rna, metadata = mdata)
DefaultAssay(seurat.rna) = 'RNA'
Idents(seurat.rna) = seurat.rna$cell_state1
degs <- FindAllMarkers(seurat.rna, assay = 'RNA', max.cells.per.ident = 500, 
                       logfc.threshold = 0.1,
                       only.pos = F, min.cells.feature = 20)
degs = data.table(degs)
degs = degs[!grepl(gene, pattern = '^MT-')]
saveRDS(degs, file = 'data/intermediate/rna/wilcox_degs_seurat_rna_Neuroglia_revision.rds')

seurat.atac = readRDS('SeuratObjects/seurat_atac_tumor_signac_RemovedCls2_4.rds')
meta_score = readRDS('MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4_withModuleScores.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = meta_score)

degs[, 'pct.diff' := pct.1 - pct.2]
degs = degs[avg_log2FC > 0.5 & pct.diff > 0]
DefaultAssay(seurat.atac) <- 'ACTIVITY'
for(state0 in unique(degs$cluster)){
  p.genes = degs[cluster == state0]$gene
  if(length(p.genes) > 50) p.genes = p.genes[1:50]
  fname = gsub('/', '.', paste0(state0, '.wilcoxDEG.Score1'), fixed = T)
  fname = gsub(' ', '.', fname, fixed = T)
  fname = gsub('+', '_', fname, fixed = T)
  seurat.atac <- AddModuleScore(seurat.atac, features = list(p.genes), 
                                name =  gsub('Score1', 'Score', fname),
                                assay = 'ACTIVITY')
  
  #fname = gsub('-', '.', fname, fixed = T)
  p1 <- FeaturePlot(seurat.atac, features = fname, max.cutoff = 'q99') +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  ggsave(p1, filename = paste0('Figures/ATAC/States/Revision/Tumor_', fname, '.pdf'),
         width = 8, height = 6, device = 'pdf')
}

saveRDS(seurat.atac@meta.data, file = 'MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4_withModuleScores.rds')

sele.genes = degs[cluster == 'Transition 1']$gene

seurat.rna <- ScaleData(seurat.rna, features = sele.genes)
seurat.tmp <- subset(seurat.rna, down = 100)
p0 = DoHeatmap(seurat.tmp, features = sele.genes)
ggsave(p0, file = 'Figures/ATAC/States/Revision/transition1_wilcox_degs_rna_heatmap.pdf', 
       width = 8, height = 6, device = 'pdf')


## diff TF activity (updated chromvar obj) ####
seurat.atac = readRDS('SeuratObjects/seurat_atac_tumor_signac_RemovedCls2_4.rds')
mdata = readRDS('MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4_withModuleScores.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = mdata)
seurat.atac = subset(seurat.atac, prediction.state.score.max1 > 0.4)
mdata = seurat.atac@meta.data

chromvar.obj = readRDS('chromVARObjects/chromvar_seurat_atac_tumor_signac_RemovedCls2_4.rds')
zscore.all = chromvar.obj@assays@data$z
dscore.all = chromvar.obj@assays@data$deviations
rownames(zscore.all) = rownames(dscore.all) = chromvar.obj@elementMetadata$name
shared.bc = intersect(rownames(mdata), colnames(dscore.all))
dscore.all = dscore.all[, shared.bc]
zscore.all = zscore.all[, shared.bc]

clusters = data.table('cluster' = mdata$predicted.state1, 
                      'barcode' = rownames(mdata)) ## change this if you will have different comparisons

clusters = clusters[barcode %in% shared.bc]

diff.tfs = runDiffMotifEnrich(dscore.all, clusters, topn = 30, 
                              max_cell_per_clust = 200, fdr = 0.1,
                              min_frac_per_cluster = 0.1, order_by = 'pvalue')
diff.tfs[, 'delta' := mean1 - mean0]
diff.tfs = diff.tfs[!grepl(feature, pattern = '^ENSG')]
diff.tfs = diff.tfs[order(-delta)]

saveRDS(diff.tfs, file = 'data/intermediate/enriched_TF_byTumorStates1.rds')

# for visualization
diff.tfs$cluster1 = factor(diff.tfs$cluster1, levels = c( 'AC-like', 'GPC-like',
                                                            'MES-like', 'Transition 2', 
                                                            'Cycling', 'HSP+ Cells',
                                                            'OPC/NPC-like', 'OC-like',
                                                            'NEU-like', 'Transition 1'))
diff.tfs = diff.tfs[order((diff.tfs$cluster1)), ]

#setkey(diff.tfs, cluster1)

cell_state_colors = print(scales::hue_pal()(10))
names(cell_state_colors) = sort(unique(diff.tfs$cluster1))

set.seed(2024)
sele.ids = sample(1:nrow(clusters), 1000)
col_ann = data.frame('cluster' = clusters[sele.ids, ]$cluster)
rownames(col_ann) = clusters[sele.ids, ]$barcode
col_ann$cluster = factor(col_ann$cluster, levels = c('GPC-like','Cycling',
                                                     'Transition 1',   'Transition 2',
                                                     'HSP+ Cells', 'MES-like', 
                                                     'OPC/NPC-like', 'AC-like', 'OC-like',
                                                     'NEU-like'))

col_ann = col_ann[order(col_ann$cluster), , drop = F]

sele.tfs = NULL
for(cl0 in unique(diff.tfs$cluster1)){
  tmp = diff.tfs[cluster1 == cl0]$feature
  if(length(tmp) == 0) next
  sele.tfs = unique(c(sele.tfs, tmp[1:min(10, length(tmp))]))
}

zscore.sele = zscore.all[unique(diff.tfs$feature), rownames(col_ann)]

zscore.sele[zscore.sele > quantile(zscore.sele, 0.95)] = quantile(zscore.sele, 0.95)
zscore.sele[zscore.sele < quantile(zscore.sele, 0.1)] = quantile(zscore.sele, 0.1)

#library(JLutils) # to select bet cutoff for cutree
#best.cutree(p1$tree_row, min = 3, max = 20, loss = FALSE, graph = FALSE)
ann_color = list('cluster' = c(`MES-like` = '#619CFF',  `Transition 2` = '#00BA38',
                               `AC-like` = '#F8766D', `GPC-like` = '#D39200',
                               `NEU-like` = '#DB72FB', `OPC/NPC-like` = '#FF61CC',
                               `Transition 1` = '#00C19F', `OC-like` = '#00B9E3',
                               `Cycling` = '#cccccc', `HSP+ Cells` = 'red'))
p3 <- pheatmap(zscore.sele,
               cluster_rows = T, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann, 
               annotation_colors = ann_color,
               gaps_col = cumsum(as.numeric(table(col_ann$cluster))))

ggsave(p3, filename = 'Figures/ATAC/States/enriched_TF_byState1.pdf',
       width = 10, height = 14, device = 'pdf')

##record the order of tfs
write(p3$tree_row$labels[p3$tree_row$order], file = 'Figures/ATAC/States/ordered_tfs_enriched_TF_byState1.txt')


write.table(diff.tfs, file = 'Tables/enriched_TF_byState1.tsv', sep = '\t',
            row.names = F, quote = F)

##save average tf chromvar zscore
avg.zscore.state <- sapply(unique(clusters$cluster), function(x){
  mtx0 = zscore.all[, clusters$cluster == x]
  return(apply(mtx0, 1, function(x) mean(x, trim = 0.025)))
})
avg.zscore.state <- avg.zscore.state[!grepl(rownames(avg.zscore.state), pattern = '^ENSG'), ]

saveRDS(avg.zscore.state, file = 'data/intermediate/avg_TF_zscore_byTumorStates1.rds')





## fc of peaks in tumor atac ####
DefaultAssay(seurat.atac) <- 'ATAC'
mtx = seurat.atac@assays$ATAC@counts
mtx_by_cls <- sapply(sort(unique(seurat.atac$predicted.state1)), function(x) {
  cl_data <- mtx[, seurat.atac$predicted.state1 == x]
  return(Matrix::rowSums(cl_data))
})
colnames(mtx_by_cls) = paste0('cluster_', sort(unique(seurat.atac$predicted.state1)))

mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = F, prior.count = 0)

mtx.fc <- sapply(1:ncol(mtx_by_cls.norm), function(x){
  fc = (mtx_by_cls.norm[, x] + 0.5) / (rowMeans(mtx_by_cls.norm[, -x]) + 0.5)
  return(log2(fc))
})
colnames(mtx.fc) = colnames(mtx_by_cls)

topn = 10000
sele.peaks <- lapply(1:ncol(mtx.fc), function(x) {
  sorted.fc <- sort(mtx.fc[, x], decreasing = T)
  filter.overall = names(which(mtx_by_cls.norm[, x] > 0.5))
  filter.fc = names(which(mtx.fc[, x] > 0))
  left.pks= intersect(filter.overall, filter.fc)
  top.peaks = intersect(names(sorted.fc[1:topn]), left.pks)
  top.peaks.dt = tidyr::separate(data.table('x' = top.peaks, 
                                            'log2FC' = sorted.fc[top.peaks],
                                            'cpm' = mtx_by_cls.norm[top.peaks, x]),  
                                 col = 'x',
                                 into = c('chr', 'start', 'end'))
  top.peaks.dt$cluster = gsub('cluster_', '', colnames(mtx.fc)[x])
  top.peaks.dt$start = as.numeric(top.peaks.dt$start)
  top.peaks.dt$end = as.numeric(top.peaks.dt$end)
  return(top.peaks.dt)
})

sapply(sele.peaks, nrow)

diff.peaks = do.call('rbind', sele.peaks)

diff.peaks$peak = paste(diff.peaks$chr, diff.peaks$start, diff.peaks$end, sep = '-')
diff.peaks = data.table(diff.peaks)
saveRDS(diff.peaks, file = 'data/intermediate/atac_pseudo_bulk_byTumorState1_topFC.rds')



## diff TFs chromVAR score of multiome data ####
chromvar.obj <- readRDS('chromVARObjects/chromvar_Tumor_Multiome_Only_Seurat_Annotated.rds')
multiome.obj <- readRDS('/mnt/isilon/tan_lab/sussmanj/pHGG/Multiome/Tumor_Multiome_Only_Seurat_Annotated.rds')

zscore.all = chromvar.obj@assays@data$z
dscore.all = chromvar.obj@assays@data$deviations
rownames(zscore.all) = rownames(dscore.all) = chromvar.obj@elementMetadata$name
shared.bc = intersect(colnames(multiome.obj), colnames(dscore.all))
dscore.all = dscore.all[, shared.bc]
zscore.all = zscore.all[, shared.bc]

clusters = data.table('cluster' = multiome.obj$cell_state1_redo, 
                      'barcode' = colnames(multiome.obj)) ## change this if you will have different comparisons

clusters = clusters[barcode %in% shared.bc]

diff.tfs = runDiffMotifEnrich(dscore.all, clusters, topn = 30, 
                              max_cell_per_clust = 500, fdr = 0.1,
                              min_frac_per_cluster = 0.1, order_by = 'pvalue')
diff.tfs[, 'delta' := mean1 - mean0]
diff.tfs = diff.tfs[!grepl(feature, pattern = '^ENSG')]
diff.tfs = diff.tfs[order(-delta)]

saveRDS(diff.tfs, file = 'data/intermediate/enriched_TF_multiome_byTumorStates1.rds')

# for visualization
diff.tfs$cluster1 = factor(diff.tfs$cluster1, levels = c( 'AC-like', 'GPC-like',
                                                          'MES-like', 'Transition 2', 
                                                          'Cycling', 'HSP+ Cells',
                                                          'OPC/NPC-like', 'OC-like',
                                                          'NEU-like', 'Transition 1'))
diff.tfs = diff.tfs[order((diff.tfs$cluster1)), ]

#setkey(diff.tfs, cluster1)

cell_state_colors = print(scales::hue_pal()(10))
names(cell_state_colors) = sort(unique(diff.tfs$cluster1))

set.seed(2024)
sele.ids = sample(1:nrow(clusters), 1000)
col_ann = data.frame('cluster' = clusters[sele.ids, ]$cluster)
rownames(col_ann) = clusters[sele.ids, ]$barcode
col_ann$cluster = factor(col_ann$cluster, levels = c('GPC-like','Cycling',
                                                     'Transition 1',   'Transition 2',
                                                     'HSP+ Cells', 'MES-like', 
                                                     'OPC/NPC-like', 'AC-like', 'OC-like',
                                                     'NEU-like'))

col_ann = col_ann[order(col_ann$cluster), , drop = F]

sele.tfs = NULL
for(cl0 in unique(diff.tfs$cluster1)){
  tmp = diff.tfs[cluster1 == cl0]$feature
  if(length(tmp) == 0) next
  sele.tfs = unique(c(sele.tfs, tmp[1:min(10, length(tmp))]))
}

zscore.sele = zscore.all[unique(diff.tfs$feature), rownames(col_ann)]

zscore.sele[zscore.sele > quantile(zscore.sele, 0.95)] = quantile(zscore.sele, 0.95)
zscore.sele[zscore.sele < quantile(zscore.sele, 0.1)] = quantile(zscore.sele, 0.1)

#library(JLutils) # to select bet cutoff for cutree
#best.cutree(p1$tree_row, min = 3, max = 20, loss = FALSE, graph = FALSE)
tumor_colors <- paletteer::paletteer_d("ggthemes::Classic_10", n = 10)
names(tumor_colors) = c('GPC-like','Cycling',
                         'Transition 1',   'Transition 2',
                         'HSP+ Cells', 'MES-like', 
                         'OPC/NPC-like', 'AC-like', 'OC-like',
                         'NEU-like')

ann_color = list('cluster' = c(`MES-like` = '#8C564BFF',  `Transition 2` = '#D62728FF',
                               `AC-like` = '#7F7F7FFF', `GPC-like` = '#1F77B4FF',
                               `NEU-like` = '#17BECFFF', `OPC/NPC-like` = '#E377C2FF',
                               `Transition 1` = '#2CA02CFF', `OC-like` = '#BCBD22FF',
                               `Cycling` = '#FF7F0EFF', `HSP+ Cells` = '#9467BDFF'))


p3 <- pheatmap(zscore.sele,
               cluster_rows = T, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann, 
               annotation_colors = ann_color,
               gaps_col = cumsum(as.numeric(table(col_ann$cluster))))

ggsave(p3, filename = 'Figures/ATAC/States/enriched_TF_multiome_byState1.pdf',
       width = 10, height = 14, device = 'pdf')

write.table(diff.tfs, file = 'Tables/enriched_TF_multiome_byState1.tsv', sep = '\t',
            row.names = F, quote = F)

## overlap with enriched TFs from snATAC-seq 
diff.tfs0 = fread('Tables/enriched_TF_byState1.tsv')

for(state0 in c('GPC-like','Cycling',
                'Transition 1',   'Transition 2',
                'HSP+ Cells', 'MES-like', 
                'OPC/NPC-like', 'AC-like', 'OC-like',
                'NEU-like')){
  message(paste0('Overlaps in ', state0, ' :'))
  message(paste(intersect(diff.tfs0[cluster1 == state0]$feature, 
                          diff.tfs[cluster1 == state0]$feature), collapse = ','))
}

## plot heatmap with top enriched TFs from snATAC

ordered.tfs = fread('Figures/ATAC/States/ordered_tfs_enriched_TF_byState1.txt', header = F)
zscore.sele = zscore.all[unique(union(diff.tfs$feature, ordered.tfs$V1)), rownames(col_ann)]

zscore.sele[zscore.sele > quantile(zscore.sele, 0.95)] = quantile(zscore.sele, 0.95)
zscore.sele[zscore.sele < quantile(zscore.sele, 0.1)] = quantile(zscore.sele, 0.1)

zscore.sele = zscore.sele[ordered.tfs$V1, ]

p4 <- pheatmap(zscore.sele,
               cluster_rows = F, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann, 
               annotation_colors = ann_color,
               gaps_col = cumsum(as.numeric(table(col_ann$cluster))))
ggsave(p4, filename = 'Figures/ATAC/States/enriched_TFfromsnATAC_multiome_byState1.pdf',
       width = 10, height = 14, device = 'pdf')


