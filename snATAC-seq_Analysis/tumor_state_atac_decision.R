source('scDataAnalysis_Utilities_simp.R')
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


## explore label transfer result ####
seurat.atac = readRDS('SeuratObjects/seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
mdata <- readRDS('MetaData/metadata_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = mdata)
Idents(seurat.atac) = seurat.atac$ATAC_snn_res.0.75

seurat.atac0 = subset(seurat.atac, prediction.state.score.max > 0.6)

cluster_label <- table(seurat.atac0$predicted.state, seurat.atac0$seurat_clusters)
cluster_label_norm <- pearson_residual(cluster_label)
pvs = pnorm(cluster_label_norm, lower.tail = F)
plot_mat <- -log10(pvs)
plot_mat[is.na(plot_mat)] = 0
plot_mat[is.infinite(plot_mat)] = max(plot_mat[!is.infinite(plot_mat)]) + 2
plot_mat_norm <- plot_mat %*% diag(1/colSums(plot_mat + 0.01)) 
colnames(plot_mat_norm) = colnames(plot_mat)

pheatmap(plot_mat_norm)

cluster_label_norm[cluster_label_norm > 50] = 50
p0 <- pheatmap(cluster_label_norm)

aa = table(seurat.atac0$seurat_clusters)/table(seurat.atac$seurat_clusters)
barplot(sort(aa), col = 'lightblue')
abline(h = 0.05, lty = 2)

ggsave(p0, filename = 'Figures/ATAC/States/matching_atac_by_cluster_FromLabelTransfer_PearsonResidual_matrix.eps',
       device = 'eps', height = 6, width = 7)

## plot signature scores ####
### < using single cell DEG genes -- not integrated ####
degs = readRDS(file = 'data/intermediate/rna/degs_seurat_rna_Neuroglia_final_noIntegration_regrMITO.rds')
seurat.atac = readRDS('SeuratObjects/seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')
degs = degs[avg_log2FC > 0.25]

#degs$score = degs$avg_log2FC
#degs = degs[order(-score)]

DefaultAssays(seurat.atac) <- 'ACTICITY'
for(state0 in unique(degs$cluster)){
  p.genes = degs[cluster == state0]$gene
  if(length(p.genes) > 50) p.genes = p.genes[1:50]
  fname = gsub('/', '.', paste0(state0, '.scDEG.Score1'), fixed = T)
  fname = gsub(' ', '.', fname, fixed = T)
  seurat.atac <- AddModuleScore(seurat.atac, features = list(p.genes), 
                                name =  gsub('Score1', 'Score', fname),
                                assay = 'ACTIVITY')
  
  #fname = gsub('-', '.', fname, fixed = T)
  p1 <- FeaturePlot(seurat.atac, features = fname, max.cutoff = 'q99') +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  ggsave(p1, filename = paste0('Figures/ATAC/States/highResl/Neuroglia_', fname, '.pdf'),
         width = 8, height = 6, device = 'pdf')
}
mdata = subset(seurat.atac@meta.data, select = c('NEU-like.scDEG.Score1', 
               'OPC.NPC-like.scDEG.Score1', 'AC-like.1.scDEG.Score1', 'AC-like.2.scDEG.Score1',
               'MES-like.scDEG.Score1', 'Interm.1.scDEG.Score1', 'Interm.2.scDEG.Score1',
               'Interm.3.scDEG.Score1', 'Cycling.scDEG.Score1'))

saveRDS(mdata, file = 'MetaData/metadata_seurat_atac_Neuroglia_signature_moduleScores.rds')

## plot signature score using all bulk limma-voom degs ####
degs = readRDS(file = 'data/intermediate/rna_pseudo_bulk_limma_voom_filterMore.rds')
degs = do.call('rbind', degs)
degs$cluster = as.character(degs$cluster)
degs$cluster = gsub('state', '', degs$cluster)
seurat.atac = readRDS('SeuratObjects/seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds')

for(state0 in unique(degs$cluster)){
  p.genes = degs[cluster == state0]$gene
  if(length(p.genes) > 50) p.genes = p.genes[1:50]
  seurat.atac <- AddModuleScore(seurat.atac, features = list(p.genes), 
                                name =  paste0(state0, '.lv.Score'),
                                assay = 'ACTIVITY')
  fname = gsub('/', '.', paste0(state0, '.lv.Score1'), fixed = T)
  fname = gsub(' ', '.', fname, fixed = T)
  fname = gsub('-', '.', fname, fixed = T)
  p1 <- FeaturePlot(seurat.atac, features = fname, max.cutoff = 'q99') +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  ggsave(p1, filename = paste0('Figures/ATAC/States/highResl/Neuroglia_', fname, '.pdf'),
         width = 8, height = 6, device = 'pdf')
}

#saveRDS(seurat.atac@meta.data, file = 'MetaData/metadata_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None_withHighReslCellStateScore.rds')

## plot signature score using all bulk shared degs and annotate atac cell states ####
degs = readRDS(file = 'data/intermediate/degs_shared_top100_filterMore.rds')

for(state0 in unique(degs$cluster)){
  p.genes = degs[cluster == state0]$gene
  fname = gsub('/', '.', paste0(state0, '.Score1'), fixed = T)
  fname = gsub(' ', '.', fname, fixed = T)
  seurat.atac <- AddModuleScore(seurat.atac, features = list(p.genes), 
                                name = gsub('Score1', 'Score', fname),
                                assay = 'ACTIVITY')
 
  p1 <- FeaturePlot(seurat.atac, features = fname, max.cutoff = 'q99') +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  ggsave(p1, filename = paste0('Figures/ATAC/States/highResl/Neuroglia_', fname, '.pdf'),
         width = 8, height = 6, device = 'pdf')
}


## determine final atac cell state ####
seurat.atac$cell_state = ''
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('9')] = 'AC-like 1'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('2', '0')] = 'AC-like 2'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('5', '7', '10')] = 'Interm 2'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('12')] = 'MES-like'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('11', '13')] = 'NEU-like'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('3', '6')] = 'OPC/NPC-like'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('1')] = 'Interm 1'  
seurat.atac$cell_state[seurat.atac$ATAC_snn_res.0.75 %in% c('4', '8', '14')] = 'Interm 3'

p2 <- DimPlot(seurat.atac, group.by = 'cell_state', label = T, label.size = 5)
saveRDS(seurat.atac@meta.data, file = "MetaData/metadata_seurat_atac_Neuroglia_signacIntegrated_bySample_withReference_None.rds")
ggsave(p2, filename = 'Figures/ATAC/States/umap_Neuroglia_cell_state_final.pdf',
       width = 8, height = 6, device = 'pdf')
