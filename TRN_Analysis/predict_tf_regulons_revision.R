## predict tf regulon per tumor state

library(data.table)
library(Seurat)
library(magrittr)
library(ggplot2)
library(matrixStats)
library(viridis)
library(openxlsx)
`%notin%` = Negate(`%in%`)

## prepare and load inputs ####

## load data for further filtering
diff.tfs = readRDS('data/intermediate/enriched_TF_byTumorStates1.rds')
daps = readRDS('data/intermediate/atac_pseudo_bulk_byTumorState_limma_voom.rds')
fc_peaks = readRDS('data/intermediate/atac_pseudo_bulk_byTumorState1_topFC.rds')
degs = readRDS('data/intermediate/rna/wilcox_degs_seurat_rna_Neuroglia_revision.rds')
degs_time = readRDS('data/intermediate/rna/degs_seurat_rna_tumor_btwTimepoint_byState1.rds')
degs_bulk = readRDS('data/intermediate/rna_pseudo_bulk_limma_voom_filterMore.rds')
motif_pk_match = readRDS('MetaData/peak_motif_match_pvE-06.rds')
motif_pk_match = as.matrix(motif_pk_match)
predicted_links = fread('EP_Prediction/regrRes4_gene_peak_links_metacell.tsv')


## prepare & filter weak peaks
if(F){
  ## just run once
  seurat.atac = readRDS('SeuratObjects/seurat_atac_tumor_signac_RemovedCls2_4.rds')
  mdata <- readRDS('MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4_withModuleScores.rds')
  seurat.atac = AddMetaData(seurat.atac, metadata = mdata)
  seurat.atac = subset(seurat.atac, prediction.state.score.max1 > 0.4)
  mtx = seurat.atac@assays$ATAC@counts
  mtx_bulk = sapply(unique(seurat.atac$predicted.state1), function(x) {
    return(rowMeans(mtx[, seurat.atac$predicted.state1==x] > 0))
  })
  pm = rowMaxs(mtx_bulk)
  weak.pks = rownames(mtx)[pm < 0.04]
  saveRDS(weak.pks, file = 'MetaData/weak_peaks_revision.rds')
  
  #peak accessible frac in cell state
  saveRDS(mtx_bulk[pm >= 0.04, ], file = 'MetaData/peak_openFrac_byState_revision.rds')
}
weak.pks = readRDS('MetaData/weak_peaks_revision.rds')
open.frac.pks = readRDS('MetaData/peak_openFrac_byState_revision.rds')

## prepare & filter unexpressed tfs
if(F){
  ## just run once
  seurat.rna = readRDS('SeuratObjects/seurat_rna_Neuroglia_final.rds')
  mdata = readRDS('MetaData/Metadata_cell_state1.rds')
  seurat.rna = AddMetaData(seurat.rna, metadata = mdata)
  
  mtx = seurat.rna@assays$RNA@counts
  expr.frac.state = sapply(unique(seurat.rna$cell_state1), function(x){
    return(rowMeans(mtx[, seurat.rna$cell_state1 == x] > 0))
  })
  colnames(expr.frac.state) = unique(seurat.rna$cell_state1)
  saveRDS(expr.frac.state, file = 'MetaData/expr_frac_state1.rds')
  
  mtx = seurat.rna@assays$RNA@data
  expr.avg.state = sapply(unique(seurat.rna$cell_state1), function(x){
    return(rowMeans(mtx[, seurat.rna$cell_state1 == x]))
  })
  colnames(expr.avg.state) = unique(seurat.rna$cell_state1)
  saveRDS(expr.avg.state, file = 'MetaData/avg_expr_state1.rds')
  
  mtx = seurat.rna@assays$RNA@data
  expr.median.state = sapply(unique(seurat.rna$cell_state1), function(x){
    return(sparseMatrixStats::rowMedians(mtx[, seurat.rna$cell_state1 == x]))
  })
  rownames(expr.median.state) = rownames(mtx)
  colnames(expr.median.state) = unique(seurat.rna$cell_state1)
  
  saveRDS(expr.median.state, file = 'MetaData/median_expr_state1.rds')
  
  ## diff expressed TFs
  avg.zscore.state = readRDS('data/intermediate/avg_TF_zscore_byTumorStates1.rds')
  tfs = rownames(avg.zscore.state)
  tfs.deg = degs[gene %in% tfs & avg_log2FC > 0.5 & p_val_adj < 0.05 & pct.1 > pct.2]
  saveRDS(tfs.deg, file = 'MetaData/tfs_deg_state1.rds')
}
expr.frac.state = readRDS('MetaData/expr_frac_state1.rds')
expr.median.state = readRDS('MetaData/median_expr_state1.rds')
expr.avg.state = readRDS('MetaData/avg_expr_state1.rds')
tfs.deg = readRDS('MetaData/tfs_deg_state1.rds')

## controlling parameters
filtered_links = NULL
filterWeakPeak = TRUE

useDAP = FALSE ## filter peaks by daps/foldChange or by peak open fraction(>0.05)
useAllDEGs = TRUE ## use all DEGs or filter
useEnhancerPeakOnly = TRUE ##use promoter enhancer interaction only or include PPIs
peakFracCutoff = 0.05 ## used to filter peaks per state (default 0.1)
tfexprfracCutoff = 0.2
if(useEnhancerPeakOnly){
  predicted_links = fread('EP_Prediction/regrRes4_gene_peak_ep_links_metacell.tsv')
}


## construct TRN per state ####
## provide or read your interested gene and tf lists
edges_list = nodes_list = list()
for(cell_state0 in colnames(expr.frac.state)){
  degs0 = degs[cluster == cell_state0 & p_val_adj < 0.01 & avg_log2FC > 0.5]
  open.frac.pks0 = open.frac.pks[, cell_state0]
  
  if(cell_state0 != 'Transition 1') degs0 = degs0[pct.1 > pct.2]
  interested_genes = unique(degs0$gene)
  
  ## as TF strength
  avg.zscore.state = readRDS('data/intermediate/avg_TF_zscore_byTumorStates1.rds')
  avg.zscore.state0 = avg.zscore.state[, cell_state0]
  
  interested_tfs = unique(c(diff.tfs[cluster1 == cell_state0]$feature, 
                              tfs.deg[cluster == cell_state0]$gene))

  ## filtering peaks
  filtered.peaks = unique(c(daps[cluster == cell_state0]$peak, 
                            fc_peaks[cluster == cell_state0]$peak))
  if(!useDAP) filtered.peaks = names(which(open.frac.pks0 > peakFracCutoff))
  
  ## further filtering TF by expression
  expr.frac.state0 <- expr.frac.state[, cell_state0]
  expr.median.state0 <- expr.median.state[, cell_state0]
  expr.avg.state0 <- expr.avg.state[, cell_state0]
  
  tf.expr.frac = expr.frac.state0[interested_tfs]
  tf.expr.frac = tf.expr.frac[!is.na(tf.expr.frac)]
  interested_tfs = names(which(tf.expr.frac > tfexprfracCutoff))

  if(length(interested_tfs) <= 5) {
    interested_tfs = intersect(names(sort(tf.expr.frac, decreasing = T))[1:10],
                               names(which(tf.expr.frac > 0.05)))
  }
  predicted_links0 = predicted_links[gene_name %in% interested_genes]
  sele.tf.mat = motif_pk_match[, interested_tfs, drop = F] 
  sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
  sele.peaks = intersect(sele.peaks, filtered.peaks)
  predicted_links0 = predicted_links0[peak_name %in% sele.peaks]
  
  if(filterWeakPeak) predicted_links0 = predicted_links0[peak_name %notin% weak.pks]
  
  ##save filtered links
  predicted_links0 = predicted_links0[Estimate > 0.2 & fdr < 0.001]
  predicted_links0$cell_state = cell_state0
  filtered_links = rbind(filtered_links, predicted_links0)
  
  ## split links by tf
  tf_regulons = list()
  open.frac.pks0 = open.frac.pks[, cell_state0]
  for(TF0 in interested_tfs){
    peaks0 = names(which(sele.tf.mat[, TF0] > 0))
    ep0 = predicted_links0[peak_name %in% peaks0]
    ep0 = subset(ep0, select = c('gene_name', 'Estimate',  'fdr', 'peak_name'))
    ep0$peak_strength = open.frac.pks0[ep0$peak_name]
    #ep0 = ep0[peak_strength > 0.05]
    ep0[, 'score' := -log10(fdr)]
    ep0$TF = TF0
    ep0[, 'edge_strength' := sum(Estimate * peak_strength), by = gene_name]
    ep0[, 'Estimate' := sum(Estimate), by = gene_name]
    ep0[, 'score' := mean(score), by = gene_name]
    ep0[, 'peaks' := paste0(peak_name, collapse = ','), by = gene_name]
    ep0[, 'peak_strengths' := paste0(peak_strength, collapse = ','), by = gene_name]
    ep0[, 'npeaks' := .N, by = gene_name]
    ep0 = subset(ep0, select = c('TF', 'gene_name', 'edge_strength', 'peak_strengths',
                                 'Estimate', 'score', 'peaks', 'npeaks'))
    ## node inf
    ep0$TF_strength = avg.zscore.state0[ep0$TF]
    ep0$gene_strength = expr.avg.state0[ep0$gene_name]
    
    ## add tf gene expr 
    ep0$TF_gene_expr = expr.avg.state0[ep0$TF]
    
    tf_regulons[[TF0]] = ep0[!duplicated(ep0)]
  }
  
  #save regulons
  cell_state1 = gsub('/', '-', cell_state0, fixed = T)
  cell_state1 = gsub(' ', '-', cell_state1, fixed = T)
  cell_state1 = gsub('+', '_', cell_state1, fixed = T)
  
  tf_regulons.comb = do.call('rbind', tf_regulons)
  
  tf_regulons.comb = tf_regulons.comb[order(-edge_strength)]
  message(paste0(cell_state0, ': ', nrow(tf_regulons.comb), ' regulons'))

  nodes1 = subset(tf_regulons.comb, select = c('TF', 'TF_strength', 'TF_gene_expr'))
  nodes1$TF_strength = nodes1$TF_strength/max(nodes1$TF_strength)
  nodes1$TF_gene_expr = nodes1$TF_gene_expr/max(nodes1$TF_gene_expr)
  nodes1$node_strength = 0.5 * nodes1$TF_strength + 0.5 * nodes1$TF_gene_expr
  nodes1 = subset(nodes1, select = c('TF', 'node_strength'))
  
  nodes2 = subset(tf_regulons.comb, select = c('gene_name', 'gene_strength'))
  names(nodes1) = names(nodes2) = c('gene', 'node_strength')
  nodes2$node_strength = nodes2$node_strength/max(nodes2$node_strength)
  
  nodes1$node_type = 'TF'
  nodes2$node_type = 'gene'
  nodes2 = nodes2[gene %notin% nodes1$gene]
  nodes = rbind(nodes1, nodes2) %>% unique()
  
  #add deg information btw timepoint
  degs_time0 = degs_time[cluster == 't2' & cell_state == cell_state0]
  setkey(degs_time0, gene)
  nodes[, 'avg_log2FC_time' := degs_time0[nodes$gene]$avg_log2FC]
  nodes[is.na(avg_log2FC_time)]$avg_log2FC_time = 0
  
  if(useDAP) {
    tf_regulon_fname = paste0('EP_Prediction/TRN/Revision/useDEGs/TF_regulons_', cell_state1, '_useDAPs.tsv')
    nodes_fname = paste0('EP_Prediction/TRN/Revision/useDEGs/nodes_', cell_state1, '_useDAPs.tsv')
    tf_rank_fname = paste0('EP_Prediction/TRN/Revision/useDEGs/ConnectRank_TF_regulons_', cell_state1, '_useDAPs.pdf')
  }else{
    tf_regulon_fname = paste0('EP_Prediction/TRN/Revision/useDEGs/TF_regulons_', cell_state1, '_peakFracCutoff', peakFracCutoff,  '.tsv')
    nodes_fname = paste0('EP_Prediction/TRN/Revision/useDEGs/nodes_', cell_state1, '_peakFracCutoff', peakFracCutoff, '.tsv')
    tf_rank_fname = paste0('EP_Prediction/TRN/Revision/useDEGs/ConnectRank_TF_regulons_', cell_state1, '_peakFracCutoff', peakFracCutoff, '.pdf')
    
  }
  
  if(!useEnhancerPeakOnly){
    tf_regulon_fname = gsub('.tsv', '_IncludePPI.tsv', tf_regulon_fname)
    nodes_fname = gsub('.tsv', '_IncludePPI.tsv', nodes_fname)
    tf_rank_fname = gsub('.pdf', '_IncludePPI.pdf', tf_rank_fname)
  }
  if(useAllDEGs)  {
    tf_regulon_fname = gsub('.tsv', '_useAllDEGs.tsv', tf_regulon_fname)
    nodes_fname = gsub('.tsv', '_useAllDEGs.tsv', nodes_fname)
    tf_rank_fname = gsub('.pdf', '_useAllDEGs.pdf', tf_rank_fname)
  }
  write.table(tf_regulons.comb, file = tf_regulon_fname,
              sep = '\t', row.names = F, quote = F)
  write.table(nodes, file = nodes_fname,
              sep = '\t', row.names = F, quote = F)
  edges_list[[cell_state0]] = tf_regulons.comb
  nodes_list[[cell_state0]] = nodes
  #plot tf ranking
  pdata = data.table('TF' = names(sort(table(tf_regulons.comb$TF))),
                     'nTargets' = as.numeric(sort(table(tf_regulons.comb$TF))))
  pdata$TF = factor(pdata$TF, levels = names(sort(table(tf_regulons.comb$TF))))
  p1 <- ggplot(data = pdata, aes(x = TF, y = nTargets, fill = nTargets)) + geom_bar(stat = 'identity') + coord_flip() +
    theme_classic() + scale_fill_viridis_c(option = "magma") + ggtitle(cell_state0) + xlab('')
 
  ggsave(p1, filename = tf_rank_fname, device = 'pdf', width = 5, height = 5)
}

if(useDAP) {
  filtered_ep_file = 'EP_Prediction/filtered_regrRes4_gene_peak_links_metacell_useDAPs_ep_revisioin.tsv'
}else{
  filtered_ep_file = paste0('EP_Prediction/filtered_regrRes4_gene_peak_links_metacell_peakFracCutoff',
                            peakFracCutoff, '_ep_revision.tsv')
}
if(!useEnhancerPeakOnly) {
  filtered_ep_file = gsub('ep.tsv', 'IncludePPI.tsv', filtered_ep_file, fixed = T)
}
if(useAllDEGs)  filtered_ep_file = gsub('.tsv', '_useAllDEGs.tsv', filtered_ep_file)

fwrite(filtered_links, file = filtered_ep_file, sep = '\t')

names(edges_list)[2] = names(nodes_list)[2] = 'OPC-NPC-like'
writexl::write_xlsx(edges_list,
                    path = 'EP_Prediction/TRN/Revision/useDEGs/TRN_edges.xlsx')

writexl::write_xlsx(nodes_list,
                    path = 'EP_Prediction/TRN/Revision/useDEGs/TRN_nodes.xlsx')

## check #up#/down genes ####
states = openxlsx::getSheetNames('EP_Prediction/TRN/Revision/useDEGs/TRN_nodes.xlsx')
for(state0 in states){
  dd = read.xlsx('EP_Prediction/TRN/Revision/useDEGs/TRN_nodes.xlsx', sheet = state0)
  dd = data.table(dd)
  perc = nrow(dd[avg_log2FC_time > 0])/nrow(dd)
  nperc = nrow(dd[avg_log2FC_time < 0])/nrow(dd)
  rt = perc/nperc
  message(paste0(state0, ': Sig up percent-> ', round(perc*100,2), 
                 ' Sig down percent-> ', round(nperc*100, 2),
                 ', up/down ratio-> ', round(rt, 2)))
}


## dotplot TFs for all states ####
pdata.comb = NULL
for(cell_state0 in  c('GPC-like','Cycling',
                      'Transition 2',
                      'HSP+ Cells', 'MES-like', 
                      'OPC/NPC-like', 'AC-like', 'OC-like',
                      'NEU-like')){
  cell_state1 = cell_state0
  cell_state1 = gsub(' ', '-', cell_state1)
  cell_state1 = gsub('/', '-', cell_state1, fixed = T)
  cell_state1 = gsub('+', '_', cell_state1, fixed = T)
  
  tf_regulon_file = paste0('EP_Prediction/TRN/Revision/useDEGs/TF_regulons_', 
                           cell_state1, '_peakFracCutoff0.05_useAllDEGs', '.tsv')
  if(!file.exists(tf_regulon_file)) next
  tf_regulons = fread(tf_regulon_file)
  
  pdata = subset(tf_regulons, select = c(TF, gene_name, TF_strength)) %>% unique()
  names(pdata)[3] = 'chromv_zscore'
  pdata[, 'nTarget' := .N, by = TF]
  pdata$N = length(unique(pdata$gene_name))
  pdata[, 'frac_target' := nTarget/N]
  pdata = subset(pdata, select = c(TF, frac_target, chromv_zscore, N)) %>% unique()
  pdata = pdata[chromv_zscore > 0]
  pdata = pdata[order(-frac_target)]
  pdata = pdata[1:min(15, nrow(pdata))]
  pdata$state = cell_state0
  pdata.comb = rbind(pdata.comb, pdata)
}
pdata.comb[chromv_zscore > 4]$chromv_zscore = 4

#pdata.comb = pdata.comb[order(state)]

pdata.comb$state = factor(pdata.comb$state, levels = c('GPC-like','Cycling',
                                                        'Transition 2',
                                                       'HSP+ Cells', 'MES-like', 
                                                       'OPC/NPC-like', 'AC-like', 'OC-like',
                                                       'NEU-like'))
pdata.comb$TF = factor(pdata.comb$TF, levels = rev(unique(pdata.comb$TF)))

ph <- ggplot(data = pdata.comb, aes(x = state, y = TF, color = chromv_zscore)) +
  geom_point(aes(size = frac_target)) + 
  scale_color_viridis_c(option = "C") + theme_bw() +
  ggtitle('') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="chromv.zscore"), 
         size=guide_legend(title="targets.coverage")) 

ggsave(ph, filename = 'EP_Prediction/TRN/Revision/useDEGs/dotplot_TFs_allStates.pdf',
       device = 'pdf', width = 4, height = 8)

names(pdata.comb)[4] = 'nDegTargetByState'
saveRDS(pdata.comb, file = 'MetaData/dotplot_TF_nTargets_revision.rds')


## sele sub-network for presentation ####
for(cell_state0 in c('GPC-like','Cycling',
                     'Transition 2',
                     'HSP+ Cells', 'MES-like', 
                     'OPC/NPC-like', 'AC-like', 'OC-like',
                     'NEU-like') ){
  cell_state1 = cell_state0
  cell_state1 = gsub(' ', '-', cell_state1)
  cell_state1 = gsub('/', '-', cell_state1, fixed = T)
  cell_state1 = gsub('+', '_', cell_state1, fixed = T)
  
  
  egs = fread(paste0('EP_Prediction/TRN/Revision/useDEGs/TF_regulons_', 
                     cell_state1, '_peakFracCutoff0.05_useAllDEGs.tsv'))
  nodes = fread(paste0('EP_Prediction/TRN/Revision/useDEGs/nodes_', 
                       cell_state1, '_peakFracCutoff0.05_useAllDEGs.tsv'))
  topg = degs[cluster == cell_state0 & avg_log2FC > 0.5]$gene
  
  egs[, 'tf_degree' := .N, by = TF]
  tf_rk = subset(egs, select = c('TF', 'tf_degree')) %>% unique()
  tf_rk = tf_rk[order(-tf_degree)]
  tf_rk = tf_rk[1:min(15, nrow(tf_rk))]
  ## old version for mes-like
  #topg = degs[cluster == 'MES-like & avg_log2FC > 1 & p_val_adj < 0.05]$gene
  #tf_nodes = nodes[node_type == 'TF' & gene %in% c('RUNX1', 'JUND', 'MAFB', 'FOSL2', 'JUN')]
  
  topg = topg[1:min(50, length(topg))]
  
  tf_nodes = nodes[node_type == 'TF' & node_strength > 0 ]
  down_expr_genes = degs[cluster == cell_state0 & avg_log2FC < 0]$gene
  if(cell_state0 != 'GPC-like') tf_nodes = tf_nodes[gene %notin% down_expr_genes & gene %in% tf_rk$TF]
  
  gene_nodes =  nodes[node_type != 'TF' & gene %in% topg]
  
  sele_nodes = rbind(tf_nodes, gene_nodes)
  sele_egs = egs[TF %in% tf_nodes$gene & gene_name %in% gene_nodes$gene]
  
  write.table(sele_egs, file = paste0('EP_Prediction/TRN/Revision/useDEGs/subnetworks/sele_edges_', cell_state1, '.tsv'), 
              sep = '\t', row.names = F, quote = F)
  write.table(sele_nodes, file = paste0('EP_Prediction/TRN/Revision/useDEGs/subnetworks/sele_nodes_',cell_state1, '.tsv'), 
              sep = '\t', row.names = F, quote = F)
  
}



## save diff tf ####
diff.tfs = readRDS('data/intermediate/enriched_TF_byTumorStates1.rds')
diff.tf.list = split(diff.tfs, by = 'cluster1')
names(diff.tf.list)[6] = 'OPC-NPC-like'
writexl::write_xlsx(diff.tf.list,
                    path = 'EP_Prediction/TRN/useDEGs/differential_tf_dscore_byCellState1.xlsx')

