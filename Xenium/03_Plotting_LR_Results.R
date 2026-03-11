library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table)
library(Seurat)
library(poolr)
library(viridis)

setwd("/mnt/isilon/tan_lab/sussmanj/Xenium/pHGG")

myeloid_tumor = read.csv("Empirical_shuffling_results_Myeloid_Tumor_30_Log.csv")
#myeloid_tumor = read.csv("Empirical_shuffling_results_Myeloid_Tumor_15_Log.csv")
#myeloid_tumor = read.csv("Empirical_shuffling_results_Myeloid_Tumor_50_Log.csv")

myeloid_tumor = na.omit(myeloid_tumor)

myeloid_tumor$LR_Pair = paste0(myeloid_tumor$Ligand, "_", myeloid_tumor$Receptor)
myeloid_tumor$Adj_P_Value = p.adjust(myeloid_tumor$P_Value, method = "BH")

average_df <- myeloid_tumor %>%
  group_by(LR_Pair) %>%
  summarize(Sample = "Aggregated",
            Actual_Score = mean(Actual_Score),
            Fold_Change = mean(Fold_Change), 
            Combined_Z_Score = stouffer(P_Value)$statistic,  
            Combined_P_Value = stouffer(P_Value)$p, 
            Num_Significant_Samples = sum(Adj_P_Value < 0.05))
average_df$Adj_P_Value = p.adjust(average_df$Combined_P_Value, method = "BH")
average_df$Combined_Sig <- ifelse(average_df$Adj_P_Value < 0.05, TRUE, FALSE)

filter_df = average_df[average_df$Combined_Sig == TRUE & average_df$Num_Significant_Samples >= 2, ]
filter_average_df = filter_df[!(filter_df$LR_Pair %in% c("CD38_PECAM1", "HLA-DQA1_CD4", "HLA-DRB5_CD4", "NCAM1_PTPRZ1", 
                                                 "CNTN1_NOTCH1", "CNTN1_NOTCH2", "CNTN1_PTPRZ1", 
                                                 "NRXN1_NLGN3", "EFEMP1_EGFR", "IL6_FCGR3B", "LGALS1_PTPRC", "NTRK3_PTPRD")), ]
filter_average_df = filter_average_df[, c("LR_Pair", "Sample", "Fold_Change", "Adj_P_Value")]
filter_individual_df = myeloid_tumor[myeloid_tumor$LR_Pair %in% filter_average_df$LR_Pair, ]
filter_individual_df = filter_individual_df[!(filter_individual_df$LR_Pair %in% c("CD38_PECAM1", "HLA-DQA1_CD4", "HLA-DRB5_CD4", "NCAM1_PTPRZ1", 
                                                         "CNTN1_NOTCH1", "CNTN1_NOTCH2", "CNTN1_PTPRZ1", 
                                                         "NRXN1_NLGN3", "EFEMP1_EGFR", "IL6_FCGR3B", "LGALS1_PTPRC", "NTRK3_PTPRD")), ]
filter_individual_df = filter_individual_df[,c("LR_Pair", "Sample", "Fold_Change", "Adj_P_Value")]

x = 1e-100
combined_df = rbind(filter_individual_df, filter_average_df)
combined_df$NegLogFDR = -log10(combined_df$Adj_P_Value + x)
combined_df$target = paste0("Tumor: ", combined_df$Sample)
combined_df$source = paste0("Myeloid: ", combined_df$Sample)
combined_df <- combined_df %>%
  separate(LR_Pair, into = c("ligand", "receptor"), sep = "_")
combined_df <- combined_df %>%
  arrange(ligand)
cor.test(filter_individual_df$Fold_Change, -log10(filter_individual_df$Adj_P_Value+x), method = "spearman")
combined_df = na.omit(combined_df)

combined_df$target = factor(combined_df$target, levels = rev(c("Tumor: pHGG_161", "Tumor: pHGG_6477", "Tumor: pHGG_339", "Tumor: pHGG_942", "Tumor: Aggregated")))

dotplot_simple <- function(combined_df, size_range = c(1, 8)) {
  entities <- c("ligand", "receptor")
  liana_mod <- combined_df %>%
    unite(entities, col = "interaction", sep = " -> ") %>%
    unite(c("source", "target"), col = "source_target", remove = FALSE) %>%
    mutate(interaction = factor(interaction, levels = rev(unique(interaction))),
           across(where(is.character), as.factor))
  
  ggplot(liana_mod, aes(x = interaction, y = target, colour = Fold_Change, size = NegLogFDR)) +  # x and y swapped
    geom_point() +
    scale_color_viridis(option = "viridis", direction = 1) +  
    scale_size_continuous(range = size_range) +
    scale_x_discrete(limits = rev(levels(liana_mod$interaction))) +  
    theme_bw() + 
    theme(
      panel.grid = element_blank(), 
      axis.title = element_blank(),  
      strip.text = element_blank(),  
    ) 
}

p1 = dotplot_simple(
  combined_df,
  size_range = c(0.2, 10) 
) + RotatedAxis() 
p1
ggsave(plot = p1, filename = "Figures/Dotplot_Xenium_30.pdf", width = 12, height = 2.6)

