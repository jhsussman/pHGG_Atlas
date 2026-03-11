library(Seurat)
library(Signac)
library(data.table)
library(magrittr)
library(Matrix)
library(matrixStats)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

setwd('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/')

addLinks <- function(multiome.obj, regr.res, tss.ann){
  peaks.range = StringToGRanges(regr.res$peak_name)
  regr.res$midpoint = floor(start(peaks.range)/2 + end(peaks.range)/2)
  tss.ann.dt = data.table('chr' = as.character(seqnames(tss.ann)),
                          'start' = start(tss.ann),
                          'gene_name' = tss.ann$gene_name)
  
  setkey(tss.ann.dt, gene_name)
  regr.res = regr.res[gene_name %in% tss.ann.dt$gene_name]
  peaks.range = StringToGRanges(regr.res$peak_name)
  
  regr.res[, 'tss' := tss.ann.dt[J(regr.res$gene_name)]$start]
  
  regr.res[, 'start' := ifelse(tss < midpoint, tss, midpoint)]
  regr.res[, 'end' := ifelse(tss < midpoint, midpoint, tss)]
  regr.res[, c('tss', 'midpoint') := NULL]
  regr.res$chromosome = as.character(peaks.range@seqnames)
  links.sig <- makeGRangesFromDataFrame(df = regr.res, seqnames.field = 'chromosome', 
                                        keep.extra.columns = TRUE)
  
  multiome.obj@assays$ATAC@links <- links.sig
  return(multiome.obj)
}

## using malignant cells only ####
signac.obj <- readRDS('SeuratObjects/seurat_atac_tumor_withTFMotifs.rds') ## with fragments and motif inf
mdata <- readRDS('MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4_withModuleScores.rds')
mdata = subset(mdata, prediction.state.score.max1 > 0.4)

signac.obj = subset(signac.obj, cell_bc %in% rownames(mdata))
signac.obj = AddMetaData(signac.obj, metadata = mdata)

signac.obj$cell_state = factor(as.character(signac.obj$predicted.state1),
                                 levels = c( "GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                                             "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like"))
Idents(signac.obj) = signac.obj$cell_state

## plot gene body region ####
## by cell state
for(show.gene in c('CHRM3', 'PLA2G2A', 'FGFR1', 'FGFR2', 'AURKB')){
  p0 <- CoveragePlot(
    object = signac.obj,
    region = show.gene,
    peaks = FALSE,
    extend.upstream = 2000,
    extend.downstream = 2000
  )
  ggsave(p0, filename = paste0('Figures/ATAC/CoveragePlot/CoveragePlot_',
                               show.gene, '.pdf'),
         device = 'pdf', width = 7, height = 7)
}

## by timepoint4
signac.obj$timepoint = factor(as.character(signac.obj$timepoint),
                               levels = c( 'Initial CNS Tumor',
                                           'Progressive (Non-Autopsy)', 
                                           'Recurrence',
                                           'Progressive (Autopsy)'))
Idents(signac.obj) = signac.obj$timepoint
for(show.gene in c('CHRM3', 'PLA2G2A', 'FGFR1', 'FGFR2', 'AURKB')){
  p0 <- CoveragePlot(
    object = signac.obj,
    region = show.gene,
    peaks = FALSE,
    extend.upstream = 2000,
    extend.downstream = 2000
  )
  ggsave(p0, filename = paste0('Figures/ATAC/CoveragePlot/CoveragePlot_',
                               show.gene, '_byTimepoint4.pdf'),
         device = 'pdf', width = 7, height = 7)
}

## by timepoint2
signac.obj$timepoint2 = 't1'
signac.obj$timepoint2[signac.obj$timepoint != 'Initial CNS Tumor'] = 't2'

signac.obj$timepoint2 = factor(as.character(signac.obj$timepoint2),
                               levels = c('t1', 't2'))
Idents(signac.obj) = signac.obj$timepoint2
for(show.gene in c('CHRM3', 'PLA2G2A', 'FGFR1', 'FGFR2', 'AURKB')){
  p0 <- CoveragePlot(
    object = signac.obj,
    region = show.gene,
    peaks = FALSE,
    extend.upstream = 2000,
    extend.downstream = 2000
  ) &  scale_fill_manual(values = rev(brewer.pal(n = 3, 'Set1')[-3]))
  ggsave(p0, filename = paste0('Figures/ATAC/CoveragePlot/CoveragePlot_',
                               show.gene, '_byTimepoint2.pdf'),
         device = 'pdf', width = 7, height = 7)
}

## by cell_state * timepoint2

Idents(signac.obj) = signac.obj$timepoint2
for(show.gene in c('CHRM3', 'PLA2G2A', 'FGFR1', 'FGFR2', 'AURKB')){
  p0 <- CoveragePlot(
    object = signac.obj,
    region = show.gene,
    peaks = FALSE,
    split.by = 'predicted.state1',
    extend.upstream = 2000,
    extend.downstream = 2000
  ) 
  ggsave(p0, filename = paste0('Figures/ATAC/CoveragePlot/CoveragePlot_',
                               show.gene, '_byCellStateTimepoint2.pdf'),
         device = 'pdf', width = 7, height = 7)
}


## plot link information ####
#load('MetaData/overlap_Glial_loops_noIntercept.RData')
seurat.rna = readRDS('SeuratObjects/seurat_rna_Neuroglia_final.rds')
mdata = readRDS('MetaData/Metadata_cell_state1.rds')
seurat.rna = AddMetaData(seurat.rna, metadata = mdata)

DefaultAssay(seurat.rna) = 'RNA'

annotation = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/scGEL/data/annotation_hg38.rds')
gene.ann <- Signac:::CollapseToLongestTranscript(annotation)
tss.ann <- GenomicRanges::resize(gene.ann, width = 1, fix = 'start')
weak.pks = readRDS('MetaData/weak_peaks.rds')

## all ep links
ep_links = fread(file = 'EP_Prediction/regrRes4_gene_peak_ep_links_metacell.tsv')
ep_links$score = ep_links$Estimate
ep_links = ep_links[score > 0.1 & fdr < 0.001 & peak_name %notin% weak.pks]

peaksInRegion <- function(peakNames, gene0, tss.ann, up = 50000, down = 50000){
  peaks.df = tidyr::separate(data.table('pk' = peakNames), col = 'pk', 
                             into = c('chr', 'start', 'end'),
                             remove = F)
  peaks.df$start = as.integer(peaks.df$start)
  peaks.df$end =  as.integer(peaks.df$end)
  tss0 = tss.ann[tss.ann$gene_name == gene0, ]@ranges@start
  
  peaks.df[peaks.df$start >= (tss0 -up) & peaks.df$end <= (tss0 + down), ]$pk
}

show.gene = 'FGFR1'
ep_links0 = ep_links[gene_name==show.gene]
signac.obj <- addLinks(multiome.obj = signac.obj, regr.res = ep_links0, tss.ann)

signac.obj$state_time = paste0(signac.obj$predicted.state1, '_', signac.obj$timepoint2)
signac.obj$state_time = factor(signac.obj$state_time,
                               levels = c( "GPC-like_t1", "GPC-like_t2", "Cycling_t1", "Cycling_t2",
                                           "Transition 1_t1", "Transition 1_t2", "Transition 2_t1","Transition 2_t2",
                                           "HSP+ Cells_t1", "HSP+ Cells_t2", "MES-like_t1", "MES-like_t2",
                                           "OPC/NPC-like_t1", "OPC/NPC-like_t2", "AC-like_t1", "AC-like_t2",
                                           "OC-like_t1", "OC-like_t2", "NEU-like_t1", 'NEU-like_t2'))
Idents(signac.obj) = signac.obj$state_time

sele.peaks = ep_links[gene_name==show.gene]$peak_name
sele.peaks # NULL means no EP called for this gene
ep_links[gene_name==show.gene]$ep_dist/1000

up = 40*1000 ## bp
down = 160*1000
sele.peaks = peaksInRegion(sele.peaks, show.gene, tss.ann, up, down)
#sele.peaks = sele.peaks[c(1, 4)]
ranges.highlight = NULL
if(length(sele.peaks) > 0){
  ranges.highlight = StringToGRanges(sele.peaks) + 200
  ranges.highlight$color = 'skyblue'
}

tumor_colors <- paletteer::paletteer_d("ggthemes::Classic_10", n = 10)
tumor_colors <- as.vector(rep(tumor_colors, each = 2))

names(tumor_colors) = c("GPC-like_t1", "GPC-like_t2", "Cycling_t1", "Cycling_t2",
                        "Transition 1_t1", "Transition 1_t2", "Transition 2_t1","Transition 2_t2",
                        "HSP+ Cells_t1", "HSP+ Cells_t2", "MES-like_t1", "MES-like_t2",
                        "OPC/NPC-like_t1", "OPC/NPC-like_t2", "AC-like_t1", "AC-like_t2",
                        "OC-like_t1", "OC-like_t2", "NEU-like_t1", 'NEU-like_t2')

p1 <- CoveragePlot(
  object = signac.obj,
  region = show.gene,
  region.highlight = ranges.highlight,
  peaks = FALSE,
  ymax = 200,
  extend.upstream = up,
  extend.downstream = down
) & scale_fill_manual(values = tumor_colors)

p1

## add expression
seurat.rna$timepoint2 = 't1'
seurat.rna$timepoint2[seurat.rna$timepoint != 'Initial CNS Tumor'] = 't2'
seurat.rna$timepoint2 = factor(as.character(seurat.rna$timepoint2),
                               levels = c('t2', 't1'))

seurat.rna$state_time = paste0(seurat.rna$cell_state1, '_', 
                               seurat.rna$timepoint2)
seurat.rna$state_time = factor(seurat.rna$state_time,
                               levels = rev(c( "GPC-like_t1", "GPC-like_t2", "Cycling_t1", "Cycling_t2",
                                           "Transition 1_t1", "Transition 1_t2", "Transition 2_t1","Transition 2_t2",
                                           "HSP+ Cells_t1", "HSP+ Cells_t2", "MES-like_t1", "MES-like_t2",
                                           "OPC/NPC-like_t1", "OPC/NPC-like_t2", "AC-like_t1", "AC-like_t2",
                                           "OC-like_t1", "OC-like_t2", "NEU-like_t1", 'NEU-like_t2')))
Idents(seurat.rna) = seurat.rna$state_time


p2 <- VlnPlot(seurat.rna, features = show.gene, pt.size = 0) +
  coord_flip() + ggthemes::theme_tufte() + NoLegend() + 
  xlab('') + scale_fill_manual(values = tumor_colors) +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) 

combPlot <- p1 + (p2 / grid::textGrob('') + 
            patchwork::plot_layout(heights = c(2, 1))) + 
            patchwork::plot_layout(ncol = 2, widths = c(4, 1))

combPlot

ggsave(combPlot, filename = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_',
                             show.gene, '.pdf'),
       device = 'pdf', width = 12, height = 7)

## vlnplot of the enhancer accessibility
signac.obj <- Signac::RunTFIDF(signac.obj)
signac.obj$state_time_r = factor(as.character(signac.obj$state_time),
                                levels = rev(c( "GPC-like_t1", "GPC-like_t2", "Cycling_t1", "Cycling_t2",
                                                "Transition 1_t1", "Transition 1_t2", "Transition 2_t1","Transition 2_t2",
                                                "HSP+ Cells_t1", "HSP+ Cells_t2", "MES-like_t1", "MES-like_t2",
                                                "OPC/NPC-like_t1", "OPC/NPC-like_t2", "AC-like_t1", "AC-like_t2",
                                                "OC-like_t1", "OC-like_t2", "NEU-like_t1", 'NEU-like_t2')))

VlnPlot(signac.obj, features = sele.peaks[1], pt.size = 0.1,
              group.by = 'state_time_r') +
  coord_flip() + ggthemes::theme_tufte() + NoLegend() + 
  xlab('') + scale_fill_manual(values = tumor_colors) +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) 

VlnPlot(signac.obj, features = sele.peaks[2], pt.size = 0.1,
        group.by = 'state_time_r') +
  coord_flip() + ggthemes::theme_tufte() + NoLegend() + 
  xlab('') + scale_fill_manual(values = tumor_colors) +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) 

## statistical test for RNA
pv.gene = c()
mtx = seurat.rna@assays$RNA@data
mdata = seurat.rna@meta.data
for(state0 in c( "GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                 "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")){
  mtx0 = mtx[, seurat.rna$cell_state1 == state0]
  mdata0 = mdata[seurat.rna$cell_state1 == state0, ]
  pv.gene[state0] = wilcox.test(mtx0[show.gene, mdata0$timepoint2 == 't2'], 
              mtx0[show.gene, mdata0$timepoint2 == 't1'], 
              alternative = 'greater')$p.value
}
pv.gene

#peak1
pv.pk = list()
mtx = signac.obj@assays$ATAC@data
mdata = signac.obj@meta.data
for(pk0 in sele.peaks){
  pv = c()
  for(state0 in c( "GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                   "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")){
    mtx0 = mtx[, signac.obj$predicted.state1 == state0]
    mdata0 = mdata[signac.obj$predicted.state1 == state0, ]
    pv[state0] = wilcox.test(mtx0[pk0, mdata0$timepoint2 == 't2'], 
                                   mtx0[pk0, mdata0$timepoint2 == 't1'])$p.value/2
  }
  pv.pk[[pk0]] = pv
}
pv.pk

write('Saved pvalue - RNA (t2 vs t1):\n', 
      file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))
write.table(data.frame(pv.gene),  append = T, quote = F, 
      file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))

write('Saved pvalue - ATAC (t2 vs t1):\n', append = T,
      file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))

write.table(data.frame(pv.pk),  append = T, quote = F, 
            file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))

write('Saved binding TFs - ATAC (t2 vs t1):\n', append = T,
      file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))

motif_pk_match = readRDS('MetaData/peak_motif_match_pvE-06.rds')
motif_pk_match = as.matrix(motif_pk_match)
tfs = apply(motif_pk_match[sele.peaks, , drop = F], 1, function(x) names(which(x)))

for(pk0 in names(tfs)){
  write(paste0('binding TFs @', pk0), append = T,
        file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))
  
  write.table(data.frame(tfs[[pk0]]),  append = T, quote = F, 
              file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))
  write('\n', append = T,
        file = paste0('Figures/ATAC/CoveragePlot/EP_CoveragePlot_', show.gene, '.txt'))
  
  
}
