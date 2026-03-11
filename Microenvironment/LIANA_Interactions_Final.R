library(liana)
library(CellChat)
library(Seurat)
library(stringr)
library(tidyverse)
library(entropy)
library(magrittr)
library(circlize)
library(gridExtra) 
library(ComplexHeatmap) 
library(dplyr)
setwd("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Interactions/Final_LIANA/")

#Read in Seurat object and sub-annotated subset for tumor cells and myeloid samples 
seurat.liana <- readRDS(file = "/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Final_Cohort_All_snRNA-seq.RDS")
DefaultAssay(seurat.liana) <- "RNA"
ng_rna <- readRDS("../../seurat_rna_Neuroglia_final.rds")
ng_rna <- AddMetaData(ng_rna, metadata=readRDS("../../metadata_seurat_rna_Neuroglia_final.rds"))
macrophage.rna.c.f <- readRDS(file = "../../Macrophages/Reanalysis_scRNA_macrophage_filtered2.rds")

ng_rna$mergelabel <- ng_rna$cell_state0
macrophage.rna.c.f$mergelabel <- macrophage.rna.c.f$cell_labels1
nontumor <- subset(seurat.liana, subset=merged_cellType=="Neuroglial", invert=TRUE)
nontumor <- subset(nontumor, subset=merged_cellType=="Macrophage/Microglia", invert=TRUE)
nontumor$mergelabel <- nontumor$merged_cellType
table(nontumor$merged_cellType)
DefaultAssay(nontumor) <- "RNA"
DefaultAssay(ng_rna) <- "RNA"
DefaultAssay(macrophage.rna.c.f) <- "RNA"

interesting.clusters <- merge(ng_rna, y=c(macrophage.rna.c.f, nontumor))
table(interesting.clusters$mergelabel)

Idents(interesting.clusters) <- "timepoint"
interesting.clusters <- RenameIdents(interesting.clusters, "Initial CNS Tumor" = "Initial Tumor", "Recurrence" = "Timepoint_2", "Progressive (Non-Autopsy)" = "Timepoint_2", 
                                     "Progressive (Autopsy)" = "Timepoint_2")
interesting.clusters$timemerge <- Idents(interesting.clusters)

saveRDS(interesting.clusters, "Seurat_for_LIANA.rds")

#Run LIANA#
seurat.liana <- interesting.clusters
seurat.liana <- readRDS(file = "Seurat_for_LIANA.rds")

DefaultAssay(seurat.liana) <- "RNA"
Idents(seurat.liana) <- "mergelabel"
liana.output <- liana_wrap(seurat.liana, idents_col = 'mergelabel')  
liana.output.aggregate <- liana.output %>% liana_aggregate()
saveRDS(liana.output, "liana.output.alltime.allLR.RDS")
saveRDS(liana.output.aggregate, "liana.alltime.allLR.aggregate.RDS")

#Run next on each timepoint
Idents(seurat.liana) <- "timepoint"
seurat.timepoint2 <- subset(x=seurat.liana, subset = timemerge == "Timepoint_2")
seurat.timepoint1 <- subset(x=seurat.liana, subset = timemerge == "Initial Tumor")
Idents(seurat.timepoint2) <- "mergelabel"
Idents(seurat.timepoint1) <- "mergelabel"

DefaultAssay(seurat.timepoint1) <- "RNA"
liana.output <- liana_wrap(seurat.timepoint1, idents_col = 'mergelabel')  
liana.output.aggregate <- liana.output %>% liana_aggregate()
saveRDS(liana.output, "liana.output.time1.allLR.RDS")
saveRDS(liana.output.aggregate, "liana.time1.allLR.aggregate.RDS")

DefaultAssay(seurat.timepoint2) <- "RNA"
liana.output <- liana_wrap(seurat.timepoint2, idents_col = 'mergelabel')  
liana.output.aggregate <- liana.output %>% liana_aggregate()
saveRDS(liana.output, "liana.output.time2.allLR.RDS")
saveRDS(liana.output.aggregate, "liana.time2.allLR.aggregate.RDS")


liana.output.aggregate = readRDS("liana.alltime.allLR.aggregate.RDS")
write.table(liana.output.aggregate, file = "liana.alltime.allLR.aggregate.tsv", sep = '\t', quote = F)

# Make Frequency Heatmap
library(paletteer)
library(dplyr)
tumor_colors <- paletteer_d("ggthemes::Classic_10", n = 9)
tumor_colors <- as.vector(tumor_colors)
tumor_levels <-  c("AC-like 1", "AC-like 2", "Cycling", "Interm 1", "Interm 2", "Interm 3", 
                   "MES-like", "OPC/NPC-like", "NEU-like")
names(tumor_colors) <- tumor_levels
myeloid_colors <- paletteer_d("ggsci::springfield_simpsons", n = 11)
myeloid_colors <- as.vector(myeloid_colors)
myeloid_levels <- c("Pre-Active MG", "Homeostatic MG", "Undetermined MG", "BMD TAM 1", "BMD TAM 2", 
                    "Lipid-Associated TAMs", "IFN-Responsive TAM", "Inflammatory TAM", 
                    "Pro-Angiogenic TAM", "Dendritic Cells", "Proliferating Myeloid")
names(myeloid_colors) <- myeloid_levels

liana_trunc <- liana.output.aggregate %>%
  filter(aggregate_rank <= 0.05) 

liana_trunc <- liana.output.aggregate %>%
  dplyr::filter(source %in% c(myeloid_levels, tumor_levels)) %>%
  dplyr::filter(target %in% c(myeloid_levels, tumor_levels)) %>%
  filter(aggregate_rank <= 0.05) 

liana_trunc$source <- factor(liana_trunc$source, levels = c(tumor_levels, myeloid_levels))
liana_trunc$target <- factor(liana_trunc$target, levels = c(tumor_levels, myeloid_levels))
heat_freq(liana_trunc, pallette=c("white", "red")) 


chord_freq(liana_trunc,
           source_groups = myeloid_levels,
           target_groups = tumor_levels) 
chord_freq(liana_trunc,
           source_groups = tumor_levels,
           target_groups = myeloid_levels) 
chord_freq(liana_trunc) 

ligand_complex <- c(
   "ANXA1", "ANXA1", "AREG", "AREG", "ARF1", "BCAN", "BCAN", "CALM1", "CALM1",
  "CALM1", "CALM1", "CALM1", "CALM1", "CALM3", "CALM3", "CALM3", "CXCL10", "CXCL10", 
  "CXCL16", "CXCL8", "CXCL8", "CXCL8", "DSCAM", 
  "FARP2", "FGF12", "FGF13", "FGF13", "FGL1", "FN1", "FN1", "FN1",  "HBEGF", "HBEGF", "HBEGF", "HBEGF", "HBEGF", "ICAM1", "ICAM1", "IGF1", 
  "IGF1", "IGF1", "IL18", "ITGAV", "ITGB2", "MAML2", "MAML2", 
  "MAML2", "MDK", "MMP2", "NAMPT", "NAMPT", "NLGN1", "NRG3",  "OSM", "PDGFB", 
  "PDGFB", "PDGFB", "PKM", "PTGS2", "PTN", "PTN", "PTN", "PTN", "SPP1", "SPP1", 
  "SPP1", "SPP1", "SPP1", "SPP1", "SPP1", "TGFA", "TGFA", "TGFB1", "TGFB1", "TGFB1", 
  "TGFB1", "TGFB1", "TGFB1", "THBS1", "THBS1", "THBS1", "THBS1", "THBS1", "TIMP1", 
  "TIMP1", "TIMP2", "TIMP2", "TIMP2", "TNC", "TNC", "TNF", "TNF", "TNF", "TNF", "TNF", 
  "VCAN", "VCAN", "VCAN","VEGFA", "VEGFA", "VEGFA", "VEGFA", "VEGFA", "VEGFA", 
  "VEGFA", "VIM", "WNT2B", "WNT5A", "WNT5A", "WNT5A", "WNT5A", "WNT5A"
)
uniqueL = unique(ligand_complex)

receptor_complex <- c(
  "FGFR1", "EGFR", "GRM7", "EGFR_ERBB2", "ERBB3", "CHRM3", "EGFR", "NRCAM", "CACNA1C", "GRM5",
  "KCNQ5", "PTH2R", "CACNA1C", "EGFR", "KCNQ5", "GRM7", "SDC4", "GRM7","DCC", 
  "ADORA1", "ADORA2B", "LIFR",  "FGFR1", "EGFR", "FGFR1", "EGFR", "CD44", 
  "ITGA3_ITGB1", "ITGAV_ITGB8", "ADCY1", "ADCY8", "ADCY9", "EDNRB", "EGFR", 
  "F2R", "IGF1R", "ADCY8", "ADCY9", "CD44", "CD9", "EGFR_ERBB2", "ERBB2_ERBB4", "PRLR", "CAV1", 
  "EGFR", "IGF1R", "INSR", "ITGA6_ITGB4", "IL1RAPL1", "THY1", "THY1",  "ITGA3", "CD44", 
  "CD47", "COLEC12", "LRP1", "MET", "MRC2", "PTPRK", "EGFR", "BOC", "NOTCH1", "NOTCH2", "NOTCH3", 
  "PTPRZ1", "FGFR1", "INSR", "ITGA5_ITGB1", "NRXN3", "EGFR", "ROBO2", "IL6ST_LIFR", "ADGRV1", 
  "ITGAV", "LRP1", "CD44", "CAV1", "CAV1", "ALK", "PTPRS", "PTPRZ1", "EGFR", "SEMA4D", 
  "MET", "PLXNA2", "PLXNA4", "DCC",  "FGFR1", "CD44", "ITGA4_ITGB1", "ITGA5_ITGB1", 
  "ITGA8_ITGB1", "ITGA9_ITGB1", "ITGAV_ITGB1", "ITGAV_ITGB5", "ERBB4", "PTPRD", "CAV1", "EGFR", 
  "ITGB1", "ITGB5", "ITGB8", "LPP", "SMAD3", "ITGA3_ITGB1", "LRP1", "LRP5", "PTPRJ", 
  "CD63", "FGFR2", "CD44", "ITGA3", "ITGB1", "CNTN1", "EGFR", "ITGA7", "NOTCH1", "PTPRS", 
  "TNFRSF1A", "TNFRSF21", "TRAF2", "CD44", "EGFR", "ITGB1", "ALK", "CD44", "EGFR", "EPHB2", 
  "GPC1", "GRIN2B", "ITGB1", "NRP2", "CD44",  "ANTXR1", "FZD3", "LRP6", "PTPRK", 
  "RYK"
)
uniqueR = unique(receptor_complex)

liana_trunc_filter <- liana_trunc %>%  
  dplyr::filter(ligand.complex %in% uniqueL) %>%
  dplyr::filter(receptor.complex %in% uniqueR) %>%
  arrange(ligand.complex, receptor.complex) 
liana_trunc_filter$source <- factor(liana_trunc_filter$source, levels = myeloid_levels)

p1 <- liana_trunc_filter %>%  
  liana_dotplot(source_groups = myeloid_levels,
                target_groups = tumor_levels, size_range = c(0.5, 5)) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  RotatedAxis()
ggsave("LIANA_Supp_Figure.pdf", plot = p1, device = "pdf", width=11.5, height = 12.5)
