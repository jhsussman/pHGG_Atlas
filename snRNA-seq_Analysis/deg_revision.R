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


## update degs among time point ####
seurat.rna = readRDS('SeuratObjects/seurat_rna_Neuroglia_final.rds')
mdata = readRDS('MetaData/Metadata_cell_state1.rds')
seurat.rna = AddMetaData(seurat.rna, metadata = mdata)
DefaultAssay(seurat.rna) = 'RNA'

seurat.rna$stage = 't1'
seurat.rna$stage[seurat.rna$timepoint != 'Initial CNS Tumor'] = 't2'
degs = NULL
for(cl0 in unique(seurat.rna$cell_state1)){
  seurat0 = subset(seurat.rna, cell_state1 == cl0)
  Idents(seurat0) = 'stage'
  degs0 = FindAllMarkers(seurat0,  max.cells.per.ident = 500, 
                         logfc.threshold = 0.25, only.pos = F)
  degs0 = data.table(degs0)
  degs0$cell_state = cl0
  degs = rbind(degs, degs0)
}
degs = degs[p_val < 0.001]
saveRDS(degs, file = 'data/intermediate/rna/degs_seurat_rna_tumor_btwTimepoint_byState1.rds')


