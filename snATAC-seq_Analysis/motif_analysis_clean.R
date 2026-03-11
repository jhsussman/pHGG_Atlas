library(chromVAR)
library(Seurat)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(pheatmap)

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



## virids style, downsample, match_cell, column color by sample & cluster
plot_enrich_tf2 <- function(sele.tfs, zscore.all, bc_clusters,
                            cluster.levels = NULL, sample.levels = NULL,
                            color_cluster = NULL, color_sample = NULL,
                            up.qt = 0.95, low.qt = 0.05,
                            ndownsample = 2000, 
                            match_cell = F, reScale = F, 
                            cluster.rows = F,
                            color_style = 'virid',
                            order.within.sample = FALSE,
                            ann_bar_names = c('cluster', 'sample')){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$sample))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(cluster0 in unique(bc_clusters$cluster)){
      tmp = bc_clusters[cluster == cluster0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'sample' = bc_clusters$sample,
                          'barcode' = bc_clusters$barcode,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  names(ann_column)[1:2] = ann_bar_names
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  ## order by levels
  if(!is.null(sample.levels)) ann_column = ann_column[order(factor(ann_column[, 2],
                                                                   levels = sample.levels[sample.levels %in% ann_column[, 2]])), ]
  ann_column = ann_column[order(factor(ann_column[, 1],
                                       levels = sort(unique(ann_column[, 1])))), ]
  
  if(is.null(color_cluster)){
    nc = length(unique(ann_column[[ann_bar_names[1]]]))
    getPalette = colorRampPalette(brewer.pal(7, "Set1"))
    if(nc >= 3) color_cluster = getPalette(nc)
    if(nc < 3) color_cluster = brewer.pal(3, "Set1")[1:nc]
    names(color_cluster) = sort(unique(bc_clusters$cluster))
  }
  
  if(is.null(color_sample)){
    nsample = length(unique(ann_column[[ann_bar_names[2]]]))
    getPalette = colorRampPalette(brewer.pal(7, "Dark2"))
    if(nsample < 3) color_sample = brewer.pal(3, "Dark2")[1:nsample]
    if(nsample > 3) color_sample = getPalette(nsample)
    if(!is.null(sample.levels)){
      names(color_sample) = sample.levels
    }else{
      names(color_sample) = unique(ann_column$sample)
    }
    
  }
  
  
  
  
  if(order.within.sample){
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column[, 1],
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column[, 2],
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
  }else{
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column[, 2],
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column[, 1],
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    
    
  }
  
  ann_colors = list('cluster' = color_cluster,
                    'sample' = color_sample)
  names(ann_colors) = ann_bar_names
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}

## high resl  ####
#seurat.tf <- readRDS('SeuratObjects/seurat_tf_dscore_tumor.rds')
#seurat.atac = readRDS('SeuratObjects/seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
mdata <- readRDS('MetaData/metadata_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
#seurat.atac = AddMetaData(seurat.atac, metadata = mdata)

#dscore.all = as.matrix(seurat.tf@assays$TF@counts)
#zscore.all = as.matrix(seurat.tf@assays$TF@scale.data)
chromvar.obj = readRDS('chromVARObjects/chromvar_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
zscore.all = chromvar.obj@assays@data$z
dscore.all = chromvar.obj@assays@data$deviations
rownames(zscore.all) = rownames(dscore.all) = chromvar.obj@elementMetadata$name
shared.bc = intersect(rownames(mdata), colnames(dscore.all))
dscore.all = dscore.all[, shared.bc]
zscore.all = zscore.all[, shared.bc]
mdata0 = subset(mdata, seurat_clusters != '2')
clusters = data.table('cluster' = mdata0$cell_state, 
                      'barcode' = rownames(mdata0),
                      'sample' = mdata0$sampleID) ## change this if you will have different comparisons
clusters = clusters[barcode %in% shared.bc]
diff.tfs = runDiffMotifEnrich(dscore.all, clusters, topn = 30, 
                              max_cell_per_clust = 500, fdr = 0.001,
                              min_frac_per_cluster = 0.2, order_by = 'pvalue')
diff.tfs[, 'delta' := mean1 - mean0]
diff.tfs = diff.tfs[!grepl(feature, pattern = '^ENSG')]
diff.tfs = diff.tfs[delta > 0.01]
## further filter and rank by mean1

saveRDS(diff.tfs, file = 'data/intermediate/enriched_TF_byTumorStates.rds')

# for visualization
diff.tfs0 = diff.tfs[delta  > 0.015] # you may need different cutoff
diff.tfs0$cluster1 = factor(diff.tfs0$cluster1, levels = c( 'AC-like 1', 'Interm 1',
                                                            'MES-like', 'Interm 3', 
                                                            'AC-like 2', 'Interm 2',
                                                            'OPC/NPC-like', 'NEU-like'))
diff.tfs0 = diff.tfs0[order((diff.tfs0$cluster1)), ]

#setkey(diff.tfs0, cluster1)

cell_state_colors = print(scales::hue_pal()(9))
names(cell_state_colors) = sort(c('Cycling', 'MES-like', 'Interm 3',
                                  'AC-like 1', 'Interm 1',
                                  'AC-like 2', 'Interm 2',
                                  'OPC/NPC-like', 'NEU-like'))

if(F){
  ph.py1 <- plot_enrich_tf2(unique(diff.tfs0$feature), dscore.all, 
                            bc_clusters = clusters,
                            up.qt = 0.95,
                            ndownsample = 1000, match_cell = F, cluster.rows = F,
                            color_cluster = cell_state_colors[unique(clusters$cluster)],
                            cluster.levels = (c('MES-like', 'Interm 3', 
                                                'AC-like 1', 'Interm 1',
                                                'AC-like 2', 'Interm 2',
                                                'OPC/NPC-like', 'NPC-like')), 
                            sample.levels = NULL,
                            reScale = T, color_style = 'purple-yellow',
                            ann_bar_names = c('cluster', 'sampleID'))
  
}


## replot

set.seed(2022)
sele.ids = sample(1:nrow(clusters), 1000)
col_ann = data.frame('cluster' = clusters[sele.ids, ]$cluster)
rownames(col_ann) = clusters[sele.ids, ]$barcode
col_ann$cluster = factor(col_ann$cluster, levels = c('MES-like',  'Interm 1',
                                                     'AC-like 1', 'AC-like 2', 
                                                     'NEU-like', 'OPC/NPC-like', 
                                                     'Interm 2', 'Interm 3'))

col_ann = col_ann[order(col_ann$cluster), , drop = F]

zscore.sele = zscore.all[unique(diff.tfs0$feature), rownames(col_ann)]

zscore.sele[zscore.sele > quantile(zscore.sele, 0.95)] = quantile(zscore.sele, 0.95)
zscore.sele[zscore.sele < quantile(zscore.sele, 0.1)] = quantile(zscore.sele, 0.1)
p2 <- pheatmap(zscore.sele,
               cluster_rows = F, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann)

#library(JLutils) # to select bet cutoff for cutree
#best.cutree(p1$tree_row, min = 3, max = 20, loss = FALSE, graph = FALSE)
ann_color = list('cluster' = c(`MES-like` = '#619CFF',  `Interm 1` = '#00BA38',
                               `AC-like 1` = '#F8766D', `AC-like 2` = '#D39200',
                               `NEU-like` = '#DB72FB', `OPC/NPC-like` = '#FF61CC',
                               `Interm 2` = '#00C19F', `Interm 3` = '#00B9E3'))
p3 <- pheatmap(zscore.sele,
               cluster_rows = T, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann, 
               annotation_colors = ann_color,
               gaps_col = cumsum(as.numeric(table(col_ann$cluster))))
ggsave(p3, filename = 'Figures/ATAC/States/Neuroglia_enriched_TF_byHighReslCluster.pdf',
       width = 10, height = 12, device = 'pdf')


## cluster TFs --not used
if(F){
  labels = cutree(p3$tree_row, k=6)
  row_ann = data.frame('group' = as.character(labels))
  rownames(row_ann) = names(labels)
  
  bks1 = seq(quantile(zscore.sele, 0.1), 0, length = 50)
  bks2 = seq(0, quantile(zscore.sele, 0.95), length = 50)
  
  p4 <- pheatmap(zscore.sele,
                 cluster_rows = T, cluster_cols = F, 
                 show_rownames = T, show_colnames = F,
                 angle_col = 45,
                 annotation_col = col_ann,
                 cutree_rows = 6, 
                 annotation_row = row_ann,
                 breaks = sort(unique(c(bks1, bks2))))
  
  table(sort(cutree(p4$tree_row, k=6)))
  
  
}

## plot pseudo-bulk matrix
zscore.sele = zscore.all[unique(diff.tfs0$feature), ]

zscore.sele.bulk <- sapply(c('MES-like',  'Interm 1',
                             'AC-like 1', 'AC-like 2', 
                             'NEU-like', 'OPC/NPC-like', 
                             'Interm 2', 'Interm 3'), function(x) {
                               bcs = intersect(clusters[cluster==x, ,drop = F]$barcode, shared.bc)
                               rowMeans(zscore.sele[, bcs])
                             })

colnames(zscore.sele.bulk) = c('MES-like',  'Interm 1',
                               'AC-like 1', 'AC-like 2', 
                               'NEU-like', 'OPC/NPC-like', 
                               'Interm 2', 'Interm 3')
zscore.sele.bulk0 = zscore.sele.bulk
zscore.sele.bulk0[zscore.sele.bulk > 2] = 2
zscore.sele.bulk0[zscore.sele.bulk < -2] = -2
p4 <- pheatmap(zscore.sele.bulk, angle_col = 45,
               cluster_rows = T, cluster_cols = T, show_colnames = T,
               show_rownames = T,
               breaks = sort(unique(c(seq(-2, 
                                          0, length.out = 50),
                                      seq(0, 
                                          2, length.out = 50)))),
               color = colorRampPalette(c('blue', 'white', 'yellow'))(99),)  

ggsave(p4, filename = 'Figures/ATAC/States/Neuroglia_enriched_TF_groupByState_bulk.pdf',
       width = 15, height = 15, device = 'pdf')


## TF avg.deviation.score ####
mdata <- readRDS('MetaData/metadata_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
chromvar.obj = readRDS('chromVARObjects/chromvar_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
zscore.all = chromvar.obj@assays@data$z
dscore.all = chromvar.obj@assays@data$deviations
rownames(zscore.all) = rownames(dscore.all) = chromvar.obj@elementMetadata$name
shared.bc = intersect(rownames(mdata), colnames(dscore.all))
dscore.all = dscore.all[, shared.bc]
zscore.all = zscore.all[, shared.bc]
mdata0 = subset(mdata, seurat_clusters != '2')
clusters = data.table('cluster' = mdata0$cell_state, 
                      'barcode' = rownames(mdata0),
                      'sample' = mdata0$sampleID) ## change this if you will have different comparisons
clusters = clusters[barcode %in% rownames(mdata0)]
dscore.all = dscore.all[, rownames(mdata0)]
zscore.all = zscore.all[, rownames(mdata0)]


avg.dscore.state <- sapply(unique(clusters$cluster), function(x){
  mtx0 = dscore.all[, clusters$cluster == x]
  return(apply(mtx0, 1, function(x) mean(x, trim = 0.025)))
})
avg.dscore.state <- avg.dscore.state[!grepl(rownames(avg.dscore.state), pattern = '^ENSG'), ]

avg.zscore.state <- sapply(unique(clusters$cluster), function(x){
  mtx0 = zscore.all[, clusters$cluster == x]
  return(apply(mtx0, 1, function(x) mean(x, trim = 0.025)))
})
avg.zscore.state <- avg.zscore.state[!grepl(rownames(avg.zscore.state), pattern = '^ENSG'), ]


saveRDS(avg.dscore.state, file = 'data/intermediate/avg_TF_dscore_byTumorStates.rds')
saveRDS(avg.zscore.state, file = 'data/intermediate/avg_TF_zscore_byTumorStates.rds')

## plot top TFs dscore
top.tfs0 = apply(avg.dscore.state, 2, function(x) names(sort(x, decreasing = T)[1:15]))
top.tfs0 = unique(as.vector(top.tfs0))

avg.mtx.top = avg.dscore.state[top.tfs0, ]
avg.mtx.top[avg.mtx.top > quantile(avg.mtx.top, 0.975)] = quantile(avg.mtx.top, 0.975)
avg.mtx.top[avg.mtx.top < quantile(avg.mtx.top, 0.1)] = quantile(avg.mtx.top, 0.1)
pp0 <- pheatmap(avg.mtx.top, angle_col = 45)

ggsave(pp0, filename = 'Figures/ATAC/States/Neuroglia_topTF_dscore_groupByState_bulk.pdf',
       width = 9, height = 12, device = 'pdf')

## plot top TFs zscore
#top.tfs = apply(avg.zscore.state, 2, function(x) names(sort(x, decreasing = T)[1:15]))
#top.tfs = unique(as.vector(top.tfs))

avg.mtx.top = avg.zscore.state[top.tfs0, ]
col_new_order = c('MES-like', 'Interm 1', 'AC-like 1', 'AC-like 2', 
                  'OPC/NPC-like', 'NEU-like', 'Interm 3', 'Interm 2')
pp <- pheatmap(avg.mtx.top[, col_new_order], cluster_cols = F,, angle_col = 45, 
               breaks = sort(unique(c(seq(-3, 0, length.out = 40),
                                      seq(0, 3, length.out = 50)))))
ggsave(pp, filename = 'Figures/ATAC/States/Neuroglia_topTF_zscore_groupByState_bulk.pdf',
       width = 9, height = 12, device = 'pdf')

## compare by cluster ####
seurat.atac = readRDS('SeuratObjects/seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
mdata <- readRDS('MetaData/metadata_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = mdata)
Idents(seurat.atac) = seurat.atac$ATAC_snn_res.0.75

#dscore.all = as.matrix(seurat.tf@assays$TF@counts)
#zscore.all = as.matrix(seurat.tf@assays$TF@scale.data)
chromvar.obj = readRDS('chromVARObjects/chromvar_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
zscore.all = chromvar.obj@assays@data$z
dscore.all = chromvar.obj@assays@data$deviations
rownames(zscore.all) = rownames(dscore.all) = chromvar.obj@elementMetadata$name
shared.bc = intersect(colnames(seurat.atac), colnames(dscore.all))
dscore.all = dscore.all[, shared.bc]
zscore.all = zscore.all[, shared.bc]

clusters = data.table('cluster' = seurat.atac$ATAC_snn_res.0.75, 
                      'barcode' = colnames(seurat.atac),
                      'sample' = seurat.atac$sampleID) ## change this if you will have different comparisons
clusters = clusters[barcode %in% shared.bc]
diff.tfs = runDiffMotifEnrich(dscore.all, clusters, topn = 20, 
                              max_cell_per_clust = 500,
                              min_frac_per_cluster = 0.2)
diff.tfs[, 'delta' := mean1 - mean0]
diff.tfs = diff.tfs[!grepl(feature, pattern = '^ENSG')]
diff.tfs0 = diff.tfs[delta  > 0.01] # you may need different cutoff
diff.tfs0$cluster1 = factor(diff.tfs0$cluster1, levels = paste0(0:14))
diff.tfs0 = diff.tfs0[order((diff.tfs0$cluster1)), ]

#setkey(diff.tfs0, cluster1)

cell_state_colors = print(scales::hue_pal()(15))
names(cell_state_colors) = paste0(0:14)
set.seed(2022)
sele.ids = sample(1:ncol(dscore.all), 2000)
col_ann = data.frame('cluster' = clusters[sele.ids, ]$cluster)
rownames(col_ann) = clusters[sele.ids, ]$barcode
col_ann$cluster = factor(col_ann$cluster, levels = paste0(0:14))

col_ann = col_ann[order(col_ann$cluster), , drop = F]

zscore.sele = zscore.all[unique(diff.tfs0$feature), rownames(col_ann)]

zscore.sele[zscore.sele > quantile(zscore.sele, 0.95)] = quantile(zscore.sele, 0.95)
zscore.sele[zscore.sele < quantile(zscore.sele, 0.1)] = quantile(zscore.sele, 0.1)
p2 <- pheatmap(zscore.sele,
               cluster_rows = F, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann)


p3 <- pheatmap(zscore.sele,
               cluster_rows = T, cluster_cols = F, show_colnames = F,
               show_rownames = T, annotation_col = col_ann)
ggsave(p3, filename = 'Figures/ATAC/States/Neuroglia_enriched_TF_groupByCluster.pdf',
       width = 10, height = 15, device = 'pdf')

## plot pseudo-bulk matrix
zscore.sele.bulk <- sapply(0:13, function(x) {
  bcs = rownames(col_ann[col_ann$cluster==x, ,drop = F])
  rowMeans(zscore.sele[, bcs])
})

colnames(zscore.sele.bulk) = 0:13
p4 <- pheatmap(zscore.sele.bulk, angle_col = 45,
               cluster_rows = T, cluster_cols = T, show_colnames = T,
               show_rownames = T)  

ggsave(p4, filename = 'Figures/ATAC/States/Neuroglia_enriched_TF_groupByCluster_bulk.pdf',
       width = 15, height = 15, device = 'pdf')
