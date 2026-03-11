library(openxlsx)
library(ggplot2)
library(ggridges)
`%notin%` = Negate('%in%')
library(pheatmap)

## ** explore and define tumor states using rna-seq data ** ##

seurat.rna = readRDS('SeuratObjects/seurat_rna_Neuroglia_rpca_bySample_veg2000.rds')

## Neftel celluar states gene module signatures ####
module.genes <- read.xlsx('MetaData/module_gene_lists_GBM_Suva2019Cell.xlsx', startRow = 4)

for(name0 in names(module.genes)[1:6]){
  p.genes = module.genes[, name0]
  p.genes = p.genes[!is.na(p.genes)]
  #seurat.atac <- signatureScore_zscore(seurat.atac, p.genes, vars.to.regress = 'nCount_ACTIVITY',
  #                                     score.name = paste0(name0, '.score'))
  seurat.rna <- AddModuleScore(seurat.rna, features = list(p.genes), 
                               name =  paste0(name0, '.NeftelScore'),
                                 assay = 'RNA')
  
}
seurat.rna$MES.NeftelScore1 = pmax(seurat.rna$MES1.NeftelScore1, seurat.rna$MES2.NeftelScore1)
seurat.rna$NPC.NeftelScore1 = pmax(seurat.rna$NPC1.NeftelScore1, seurat.rna$NPC2.NeftelScore1)

## plot
for(state0 in c('AC.NeftelScore1', 'MES.NeftelScore1', 
                'MES1.NeftelScore1','MES2.NeftelScore1',
                'NPC.NeftelScore1', 'OPC.NeftelScore1',
                'NPC1.NeftelScore1', 'NPC2.NeftelScore1')){
  p0 <-FeaturePlot(seurat.rna, features = state0, max.cutoff = 'q99')+ 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  ggsave(p0, filename = paste0('Figures/RNA/States/Neuroglia_Neftel_state_', state0, '.pdf'),
         device = 'pdf', width = 8, height = 6)
}



## cycle scores ####
cycle3 = fread('MetaData/regev_lab_cell_cycle_genes.txt', header = F)$V1
s.genes = cycle3[1:43]
g2m.genes = cycle3[44:97]
seurat.rna <- CellCycleScoring(seurat.rna, s.features = s.genes, 
                                 g2m.features = g2m.genes, 
                                 assay = 'RNA')

p0 <- FeaturePlot(seurat.rna, features = 'S.Score', max.cutoff = 'q99')+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p1 <- FeaturePlot(seurat.rna, features = 'G2M.Score', max.cutoff = 'q99')+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

seurat.rna$timepoint = factor(seurat.rna$timepoint, levels = 
                              c('Initial CNS Tumor', 'Progressive (Non-Autopsy)',
                                'Recurrence', 'Progressive (Autopsy)'))
seurat.rna$molecularClass = sapply(seurat.rna$molecularClass, 
                                     function(x) gsub('HGG - ', '', x))
seurat.rna$molecularClass = sapply(seurat.rna$molecularClass, 
                                     function(x) gsub(' ', '_', x))
seurat.rna$molecularClass = factor(seurat.rna$molecularClass, levels = 
                                     c('H3F3A_K27M','H3F3A_G34V',
                                       'IDH1_R132H', 'NOS'))

timeColors <- brewer.pal(4, 'Set2')
moleculeColors <- brewer.pal(5, 'Paired')[-1]

names(timeColors) = c('Initial CNS Tumor', 'Progressive (Non-Autopsy)',
                      'Recurrence', 'Progressive (Autopsy)')
names(moleculeColors) = c('H3F3A_K27M','H3F3A_G34V',
                          'IDH1_R132H', 'NOS')

p2 <- DimPlot(seurat.rna, group.by = 'timepoint') + 
  scale_color_manual(values = timeColors)
p3 <- DimPlot(seurat.rna, group.by = 'molecularClass') +
  scale_color_manual(values = moleculeColors)


ggsave(p0, filename = paste0('Figures/RNA/States/S_Neuroglia_state.pdf'),
       device = 'pdf', width = 8, height = 6)
ggsave(p1, filename = paste0('Figures/RNA/States/G2M_Neuroglia_state.pdf'),
       device = 'pdf', width = 8, height = 6)
ggsave(p2, filename = paste0('Figures/RNA/umap_Neuroglia_byTime.pdf'),
       device = 'pdf', width = 8, height = 6)
ggsave(p3, filename = paste0('Figures/RNA/umap_Neuroglia_byMolecule.pdf'),
       device = 'pdf', width = 8, height = 6)

## Wang et al 2017 cancer cell signatures ####
# getting data from sheets
sheets <- openxlsx::getSheetNames('MetaData/wang2017_CancerCell_signatures.xlsx')
module_genes <- lapply(sheets, openxlsx::read.xlsx, 
                       xlsxFile='MetaData/wang2017_CancerCell_signatures.xlsx',
                       startRow = 2)

# assigning names to data frame
names(module_genes) <- sheets

for(name0 in names(module_genes)){
  p.genes = module_genes[[name0]]$GeneSymbol
  p.genes = p.genes[!is.na(p.genes)]
  #seurat.atac <- signatureScore_zscore(seurat.atac, p.genes, vars.to.regress = 'nCount_ACTIVITY',
  #                                     score.name = paste0(name0, '.score'))
  seurat.rna <- AddModuleScore(seurat.rna, features = list(p.genes), assay = 'RNA',
                                 name =  paste0(name0, '.WangScore'))
  p0 <- FeaturePlot(seurat.rna, features = paste0(name0, '.WangScore1'),
                    max.cutoff = 'q99')+ 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  ggsave(p0, filename = paste0('Figures/RNA/States/Neuroglia_Wang_state_', name0, '.pdf'),
         device = 'pdf', width = 8, height = 6)
}

p1 <- FeaturePlot(seurat.rna, features = 'Neftel.max.score',
            max.cutoff = 'q99') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p2 <- FeaturePlot(seurat.rna, features = 'Cycle.max.score',
                  max.cutoff = 'q99') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

ggsave(p1, filename = 'Figures/RNA/States/Neuroglia_Neftel_state_max_score.pdf',
       device = 'pdf', width = 8, height = 6)
ggsave(p2, filename = 'Figures/RNA/States/Neuroglia_Cycling_max_score.pdf',
       device = 'pdf', width = 8, height = 6)
cor(seurat.rna$Neftel.max.score, seurat.rna$Cycle.max.score)


## Neftel quadrant plot ####
module.score <- data.table(subset(seurat.rna@meta.data, select = c('MES.NeftelScore1',  'AC.NeftelScore1', 
                                                                     'OPC.NeftelScore1', 'NPC.NeftelScore1', 
                                                                     'timepoint', 'molecularClass', 'neftel_type')),
                           keep.rownames = T)
module.score$y.score = pmax(module.score$NPC.NeftelScore1, module.score$OPC.NeftelScore1) - 
  pmax(module.score$MES.NeftelScore1, module.score$AC.NeftelScore1)
module.score$x.score = module.score$NPC.NeftelScore1 - module.score$OPC.NeftelScore1
module.score[y.score < 0, ]$x.score = module.score[y.score < 0, ]$MES.NeftelScore1 - 
  module.score[y.score < 0, ]$AC.NeftelScore1

module.score$Neftel_Type = 'OPC'
module.score[x.score < 0 & y.score < 0]$Neftel_Type = 'AC'
module.score[x.score > 0 & y.score < 0]$Neftel_Type = 'MES'
module.score[x.score > 0 & y.score > 0]$Neftel_Type = 'NPC'

module.score$Neftel.max.score = module.score$OPC.NeftelScore1
module.score$Neftel.max.score = ifelse(module.score$Neftel_Type == 'AC',
                                       module.score$AC.NeftelScore1, module.score$Neftel.max.score)
module.score$Neftel.max.score = ifelse(module.score$Neftel_Type == 'MES',
                                       module.score$MES.NeftelScore1, module.score$Neftel.max.score)
module.score$Neftel.max.score = ifelse(module.score$Neftel_Type == 'NPC',
                                       module.score$NPC.NeftelScore1, module.score$Neftel.max.score)


mdata = data.frame(subset(module.score, select = c('Neftel_Type', 'Neftel.max.score')))
rownames(mdata) = module.score$rn
seurat.rna <- AddMetaData(seurat.rna, metadata = mdata)
seurat.rna$Cycle.max.score = pmax(seurat.rna$G2M.Score, seurat.rna$S.Score)
DimPlot(seurat.rna, group.by = 'Neftel_Type')
DimPlot(seurat.rna, group.by = 'Phase')
FeaturePlot(seurat.rna, features = 'Neftel.max.score', max.cutoff = 'q98') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(seurat.rna, features = 'Cycle.max.score', max.cutoff = 'q98') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

saveRDS(seurat.rna, file = 'SeuratObjects/seurat_rna_Neuroglia_rpca_bySample_veg2000.rds')

## < plot by timepoint
set.seed(2022)
sele.ids = NULL
n0 = min(table(module.score$timepoint))
for(t0 in unique(module.score$timepoint)){
  ids0 = which(module.score$timepoint == t0)
  ids0 = sample(ids0, n0)
  sele.ids = c(sele.ids, ids0)
}
pdata = module.score[sele.ids, ]
pdata$timepoint = factor(pdata$timepoint, levels =  c('Initial CNS Tumor', 
                                                      'Progressive (Non-Autopsy)','Recurrence', 
                                                      'Progressive (Autopsy)'))
p0 <- ggplot(data = pdata, aes(x = x.score, y = y.score)) + 
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") + theme_classic() + facet_wrap(~ timepoint) +
  theme(legend.position = 'none') + xlab('') + ylab('')+ xlab('') + ylab('') + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 0, linetype="dashed", color = "red")
ggsave(p0, filename = 'Figures/RNA/States/rna_Neuroglia_NeftelStates_across_timepoint_quadrantPlot.eps',
       width = 6, height = 6)

## < plot by molecularClass
set.seed(2022)
sele.ids = NULL
n0 = min(table(module.score$molecularClass))
for(t0 in unique(module.score$molecularClass)){
  ids0 = which(module.score$molecularClass == t0)
  ids0 = sample(ids0, n0)
  sele.ids = c(sele.ids, ids0)
}
pdata = module.score[sele.ids, ]
p1 <- ggplot(data = pdata, aes(x = x.score, y = y.score)) + 
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") + theme_classic() + facet_wrap(~ molecularClass) +
  theme(legend.position = 'none') + xlab('') + ylab('')+ xlab('') + ylab('') + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  geom_vline(xintercept = 0, linetype="dashed", color = "red")
ggsave(p1, filename = 'Figures/RNA/States/rna_Neuroglia_NeftelStates_across_molecularClass_quadrantPlot.eps',
       width = 6, height = 6)


## 
seurat.rna <- readRDS('/mnt/isilon/tan_lab/suny6/pHGG_2022/astrocyte_integrated_final_with_ref_093022.RDS')

mdata = subset(seurat.rna@meta.data, identity == 'oligo/neural-like')

fg.freq <- table(mdata$timepoint)/nrow(mdata)
bg.freq <- table(seurat.rna$timepoint)/ncol(seurat.rna)

pdata <- data.table('freq' = c(fg.freq, bg.freq), 
                    'timepoint' = c(names(fg.freq), names(bg.freq)),
                    'condition' = c(rep('High_Invasion', 4), rep('All Tumors', 4)))


pdata$timepoint = factor(pdata$timepoint, levels = 
                                  c('Initial CNS Tumor', 'Progressive (Non-Autopsy)',
                                    'Recurrence', 'Progressive (Autopsy)'))
pdata$condition = factor(pdata$condition, levels = c('All Tumors', 'High_Invasion'))
p1 <- ggplot(data = pdata, aes(x = condition, y = freq, fill = timepoint)) +
  geom_bar(stat = 'identity', position = 'stack') + xlab('') + ylab('fraction') +
  scale_fill_manual(values = brewer.pal(4, 'Set2')) +  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(p1, filename = 'Figures/RNA/States/barplot_high_invasive_tumors_by_timepoint_rna.pdf',
       device = 'pdf', width = 4, height = 5)

## further dissect ####
seurat.rna = readRDS('/mnt/isilon/tan_lab/suny6/pHGG_2022/astrocyte_integrated_final_with_ref_093022.RDS')
DefaultAssay(seurat.rna) = 'integrated'
seurat.rna <- FindNeighbors(seurat.rna, resolution = 0.3)
DimPlot(seurat.rna, group.by = 'integrated_snn_res.0.3', label.size = 6, label = T) + NoLegend()

# reannotate cell states
seurat.rna$cell_state = 'AC-like 1'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(0, 4)] = 'AC-like 2'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(6)] = 'MES-like'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(1, 10)] = 'OPC/NPC-like'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(7)] = 'NEU-like'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(3)] = 'Interm 1'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(9)] = 'Interm 2'
seurat.rna$cell_state[seurat.rna$integrated_snn_res.0.3  %in% c(5, 8)] = 'Cycling'
p1 <- DimPlot(seurat.rna, group.by = 'cell_state', label.size = 6, label = T) 


## introduce another intermediate state
seurat.rna$cell_state0 = seurat.rna$cell_state
seurat.rna$cell_state0[seurat.rna$integrated_snn_res.0.3  %in% c(4)] = 'AC-like 2'
seurat.rna$cell_state0[seurat.rna$integrated_snn_res.0.3  %in% c(0)] = 'Interm 2'

seurat.rna$cell_state0[seurat.rna$integrated_snn_res.0.3  %in% c(3)] = 'Interm 1'
seurat.rna$cell_state0[seurat.rna$integrated_snn_res.0.3  %in% c(9)] = 'Interm 3'
p2 <- DimPlot(seurat.rna, group.by = 'cell_state0', label.size = 6, label = T) 
p2
saveRDS(seurat.rna@meta.data, file = 'MetaData/metadata_seurat_rna_Neuroglia_final.rds')
saveRDS(seurat.rna, file = 'SeuratObjects/seurat_rna_Neuroglia_final.rds')
ggsave(p2, filename = 'Figures/RNA/umap_Neuroglia_cell_state.pdf', device = 'pdf',
       width = 8, height = 6)

## update annotation for lower-resolusion version --ignore 
seurat.rna$cell_state_lowResl = seurat.rna$identity
seurat.rna$cell_state_lowResl[seurat.rna$cell_state_lowResl == 'astrocyte-reactive'] = 'AC-like 2'
seurat.rna$cell_state_lowResl[seurat.rna$cell_state_lowResl == 'astrocyte-progenitor'] = 'AC-like 1'
seurat.rna$cell_state_lowResl[seurat.rna$cell_state_lowResl == 'neural-like progenitor'] = 'NPC-like'
seurat.rna$cell_state_lowResl[seurat.rna$cell_state_lowResl == 'oligo/neural-like'] = 'OPC/NPC-like'
seurat.rna$cell_state_lowResl[seurat.rna$cell_state_lowResl == 'cycling'] = 'Cycling'
seurat.rna$cell_state_lowResl[seurat.rna$cell_state_lowResl == 'mesenchymal'] = 'MES-like'

saveRDS(seurat.rna, file = 'SeuratObjects/seurat_rna_Neuroglia_final.rds')
