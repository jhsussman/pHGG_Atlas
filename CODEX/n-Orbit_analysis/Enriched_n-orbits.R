library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

setwd("/mnt/isilon/tan_lab/yanga11/CytoIntermediates/HyperGlioma/Nolan/NolanInstances")

enrichedc1s0<-read.csv("n-orbitEnrich_Results/n-orbit-enrichmentc1s0-results.csv")
enrichedc1s0<-enrichedc1s0[,-1]
enrichedc1s2<-read.csv("n-orbitEnrich_Results/n-orbit-enrichmentc1s2-results.csv")
enrichedc1s2<-enrichedc1s2[,-1]

enrichedc1s0<-enrichedc1s0[enrichedc1s0$qvalue<0.05 & enrichedc1s0$qvalue>0,]
enrichedc1s2<-enrichedc1s2[enrichedc1s2$qvalue<0.05 & enrichedc1s2$qvalue>0,]

write.csv(enrichedc1s0, "sig-enrichc1s0.csv")
write.csv(enrichedc1s2, "sig-enrichc1s2.csv")

