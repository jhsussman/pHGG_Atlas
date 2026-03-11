library(readxl)
library(ggplot2)
library(EnhancedVolcano)
library(ggrastr)
library(S4Vectors)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)
library(dplyr)
library(gridExtra)
setwd("G:/My Drive/University of Pennsylvania/Labs/Tan Lab/HTAN/GSEA/GLM")

data <- read.csv("CPTCA_pHGG_limmavoom_mixed_model_univariate_timepoint2_allPatient.csv",
                 header = T)
dim(data)

rownames(data) = data$X
interesting_genes = c("TGFB2", "SCL12A7", "CHRM3", "RYR3", 
                      "PLA2G2A", "BCL6",  "ACHE", "CD44", "DRD2",
                      "PED10A", "OXTR", "GABRA2", "DRD2", "LIF", 
                      "MERTK", "CASP4", "KCNQ3", "NDST4", 
                      "RBP2", "CHML", "ASPN", "GNGT1", "COL2A1", 
                      "CXCL17", "WNT4", "IL7", "ARF5", "KCNQ1",
                      "HDAC10", "HDAC8", "SCN8A", "CASP1", "GABRA4", 
                      "ADORA2B", "RXRA", "IGFBP3", "MAPKAPK",
                      "NPTX2", "MOPB", "APOL3", "CNDP1", "CXCL12")

p1 <- EnhancedVolcano(data,
                      lab = rownames(data), xlab = "LogFC",
                      selectLab = interesting_genes,
                      x = 'logFC',
                      y = 'adj.P.Val',
                      pCutoff = 0.05,
                      FCcutoff = 0.5, 
                      xlim = c(-2.5, 2.6),
                      ylim = c(0, 4.6),
                      drawConnectors = TRUE,
                      labSize = 4.0, max.overlaps = 40,
                      legendPosition = "right",
                      title = "Tumor Cell GLM")
p1
p1 <- rasterize(p1, layers='Point', dpi = 600)
#ggsave("VolcanoPlot.pdf", p1, device = "pdf", width = 9, height = 7)


#GSEA

#By LogFC
res2 <- data %>% 
  dplyr::select(X, logFC) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(X) %>%
  summarize(stat=mean(logFC))

#By statistic
res2 <- data %>% 
  dplyr::select(X, t) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(X) %>%
  summarize(stat=mean(t))

res2
ranks <- deframe(res2)

pathways.go <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
pathways.kegg <- msigdbr(species = "human", category = "C2", subcategory = "KEGG")
pathways.reactome <- msigdbr(species = "human", category = "C2", subcategory = "REACTOME")
pathways.hallmark <- msigdbr(species = "human", category = "H")

pathwaysToUse <- pathways.go
pathwaysToUse <- pathways.kegg
pathwaysToUse <- pathways.reactome
pathwaysToUse <- pathways.hallmark
pathway_list_hallmark = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list_kegg = split(x = pathways.kegg$gene_symbol, f = pathways.kegg$gs_name)
pathway_list_merge = c(pathway_list_hallmark, pathway_list_kegg)

#
#neftel_modules <- read.table("Neftel_Modules.txt", header = TRUE, sep = "\t")
#neftel_list <- lapply(neftel_modules, function(x) x[x != ""])
#pathway_list_merge = neftel_list

GSEA_result <-  fgsea(pathways = pathway_list_merge, stats = ranks)  #Change this 

fgseaResTidy <- GSEA_result %>%
  as_tibble() %>%
  arrange(desc(NES))
head(GSEA_result[order(pval), ])

sig_pathways <- GSEA_result[GSEA_result$padj<0.05, ]
dim(sig_pathways)
sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
sig_pathways = na.omit(sig_pathways)

n_top <- 18
top_pathways <- head(sig_pathways, n_top)
bottom_pathways <- tail(sig_pathways, n_top)
combined_pathways <- rbind(top_pathways, bottom_pathways)
combined_pathways <- combined_pathways[order(combined_pathways$NES), ]
combined_pathways <- unique(as.data.frame(combined_pathways))
combined_pathways <- combined_pathways[!duplicated(combined_pathways[["pathway"]]), ]

combined_pathways$pathway <- factor(
  combined_pathways$pathway,
  levels = combined_pathways$pathway[order(combined_pathways$NES)]
)
p1 <- ggplot(combined_pathways, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_gradient(low = "blue", high = "yellow") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "-log10(FDR)") +  # Update axis and legend labels
  theme_minimal() + coord_flip()+
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, family = "sans", size = 10, color = "black"),
    axis.text.y = element_text(family = "sans", size = 12, color = "black"),
    axis.line = element_blank(),
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt")
  )
p1

pa = plotEnrichment(neftel_list[["MES"]],ranks) + ggtitle("MES1 + MES2") + theme_bw()
pb = plotEnrichment(neftel_list[["MES1"]],ranks) + ggtitle("MES1") + theme_bw()
pc = plotEnrichment(neftel_list[["NPC"]],ranks) + ggtitle("NPC1 + NPC2") + theme_bw()
pd = plotEnrichment(neftel_list[["NPC1"]],ranks) + ggtitle("NPC1") + theme_bw()
grid.arrange(pa, pb, pc, pd, ncol = 2)

###########MYELOID################
data <- read.csv("CPTCA_pHGG_dream_mixed_effect_model_univariate_timepoint2_myeloid_2023-10-10.csv",
                 header = T)
#GSEA
res2 <- data %>% 
  dplyr::select(X, logFC) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(X) %>%
  summarize(stat=mean(logFC))
res2
ranks <- deframe(res2)

pathways.go <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
pathways.kegg <- msigdbr(species = "human", category = "C2", subcategory = "KEGG")
pathways.reactome <- msigdbr(species = "human", category = "C2", subcategory = "REACTOME")
pathways.hallmark <- msigdbr(species = "human", category = "H")

pathwaysToUse <- pathways.go
pathwaysToUse <- pathways.kegg
pathwaysToUse <- pathways.reactome
pathwaysToUse <- pathways.hallmark
pathway_list_hallmark = split(x = pathways.hallmark$gene_symbol, f = pathways.hallmark$gs_name)
pathway_list_kegg = split(x = pathways.kegg$gene_symbol, f = pathways.kegg$gs_name)
pathway_list_merge = c(pathway_list_hallmark, pathway_list_kegg)

GSEA_result <-  fgsea(pathways = pathway_list_hallmark, stats = ranks)  #Change this 

fgseaResTidy <- GSEA_result %>%
  as_tibble() %>%
  arrange(desc(NES))
head(GSEA_result[order(pval), ])

sig_pathways <- GSEA_result[GSEA_result$padj<0.05, ]
dim(sig_pathways)
sig_pathways <- sig_pathways[order(sig_pathways$NES), ]
sig_pathways = na.omit(sig_pathways)

n_top <- 10
top_pathways <- head(sig_pathways, n_top)
bottom_pathways <- tail(sig_pathways, n_top)
combined_pathways <- rbind(top_pathways, bottom_pathways)
combined_pathways <- combined_pathways[order(combined_pathways$NES), ]
combined_pathways <- unique(as.data.frame(combined_pathways))
combined_pathways <- combined_pathways[!duplicated(combined_pathways[["pathway"]]), ]

combined_pathways$pathway <- factor(
  combined_pathways$pathway,
  levels = combined_pathways$pathway[order(combined_pathways$NES)]
)
p1 <- ggplot(combined_pathways, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_gradient(low = "blue", high = "yellow") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +  # Update axis and legend labels
  theme_minimal() + coord_flip()+
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, family = "sans", size = 10, color = "black"),
    axis.text.y = element_text(family = "sans", size = 12, color = "black"),
    axis.line = element_blank(),
    plot.margin = margin(l = 5, r = 5, b = 5, t = 5, unit = "pt")
  )
p1

