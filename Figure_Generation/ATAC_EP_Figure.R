library(Seurat)
library(Signac)
library(data.table)
library(magrittr)
library(Matrix)
library(matrixStats)
library(GenomicRanges)
library(ggplot2)
library(paletteer)

setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Analysis_Revisions")

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
signac.obj <- readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/SeuratObjects/seurat_atac_tumor_withTFMotifs.rds') ## with fragments and motif inf
mdata <- readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/MetaData/metadata_seurat_atac_tumor_signac_RemovedCls2_4_withModuleScores.rds')
mdata = subset(mdata, prediction.state.score.max1 > 0.4)

signac.obj = subset(signac.obj, cell_bc %in% rownames(mdata))
signac.obj = AddMetaData(signac.obj, metadata = mdata)

cell_states = c("GPC-like", "Cycling", "Transition 1", "Transition 2", "HSP+ Cells", 
                "MES-like", "OPC/NPC-like", "AC-like", "OC-like", "NEU-like")
signac.obj$cell_state = factor(as.character(signac.obj$predicted.state1),
                               levels = cell_states)
Idents(signac.obj) = signac.obj$cell_state

## < add link information
#load('MetaData/overlap_Glial_loops_noIntercept.RData')
seurat.rna = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/SeuratObjects/seurat_rna_Neuroglia_final.rds')
mdata = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/MetaData/Metadata_cell_state1.rds')
seurat.rna = AddMetaData(seurat.rna, metadata = mdata)

DefaultAssay(seurat.rna) = 'RNA'

annotation = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/scGEL/data/annotation_hg38.rds')
gene.ann <- Signac:::CollapseToLongestTranscript(annotation)
tss.ann <- GenomicRanges::resize(gene.ann, width = 1, fix = 'start')
weak.pks = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/MetaData/weak_peaks.rds')

## all ep links
ep_links = fread(file = '/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/EP_Prediction/regrRes4_gene_peak_ep_links_metacell.tsv')
ep_links$score = ep_links$Estimate
ep_links = ep_links[score > 0.2 & fdr < 0.001 & peak_name %notin% weak.pks]

peaksInRegion <- function(peakNames, gene0, tss.ann, up = 50000, down = 50000){
  peaks.df = tidyr::separate(data.table('pk' = peakNames), col = 'pk', 
                             into = c('chr', 'start', 'end'),
                             remove = F)
  peaks.df$start = as.integer(peaks.df$start)
  peaks.df$end =  as.integer(peaks.df$end)
  tss0 = tss.ann[tss.ann$gene_name == gene0, ]@ranges@start
  
  peaks.df[peaks.df$start >= (tss0 -up) & peaks.df$end <= (tss0 + down), ]$pk
}

show.gene = 'CHRM3'
ep_links0 = ep_links[gene_name==show.gene]
signac.obj <- addLinks(multiome.obj = signac.obj, regr.res = ep_links0, tss.ann)

sele.peaks = ep_links[gene_name==show.gene]$peak_name
sele.peaks # NULL means no EP called for this gene
tss.ann[tss.ann$gene_name == show.gene, ]

up = 100*1000 
down = 100*1000
sele.peaks = peaksInRegion(sele.peaks, show.gene, tss.ann, up, down)
ranges.highlight = NULL
if(length(sele.peaks) > 0){
  ranges.highlight = StringToGRanges(sele.peaks) + 100
  ranges.highlight$color = '#FFD300'
}

tumor_colors <- paletteer_d("ggthemes::Classic_10", n = 10)
tumor_colors <- as.vector(tumor_colors)
names(tumor_colors) = cell_states
p1 <- CoveragePlot(
  object = signac.obj,
  region = show.gene,
  region.highlight = ranges.highlight,
  peaks = FALSE,
  extend.upstream = up,
  extend.downstream = down
) & scale_fill_manual(values = tumor_colors)

## add expression

p2 <- VlnPlot(seurat.rna, features = show.gene, pt.size = 0, group.by = 'cell_state') +
  coord_flip() + ggthemes::theme_tufte() + NoLegend() + 
  xlab('') +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = tumor_colors)

combPlot <- p1 + (p2 / grid::textGrob('') + patchwork::plot_layout(heights = c(2, 1))) + 
  patchwork::plot_layout(ncol = 2, widths = c(2, 1))

combPlot


