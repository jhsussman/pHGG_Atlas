library(Clonalscope)
library(Seurat)
library(ggplot2)
library(stringr)
library(gridExtra)
library(grid)
library("cowplot")
library(dplyr)
library(readxl)
library(clevRvis)
library(fishplot)
library(RColorBrewer)

setwd("~/Desktop/Nancy_lab/CPTCA_pHGG/") 
source("analysis/plot_functions.R") # functions lineage plot
source("analysis/plot_help_functions.R") # functions for clone percentage & fishplot

# load data & metadata
snrna_obj <- readRDS("data/new/snRNA_merged_SeuratObj_QCFiltered_DoubletFiltered_downsampled_rPCA_integrated_clustered_withInferCNVTumorNormal_withNeftelTypes.rds")
metadata = read.csv("data/old/HTAN_pHGG_snRNAseq_file_annot_with_Diaz_and_Jabado_samples_2022-07-29.csv")
rownames(metadata) = metadata$fileAccession
tumor_anno=readRDS("data/new/tumor_metadata_Wenbao.rds")

# Create FishTail Plots to observe Clonal changes over time
lineage_match=read.csv("data/new/sample_match/lineage_match_progressiveSeed_14patient.txt",sep="\t",header=T)
colnames(lineage_match) = c("Patient.ID","snRNA_seed","snRNA_fileAccessions","Bulk_sample")
lineage_match = lineage_match[c(1:7,9:15),]
consistent_col = c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF","#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF")

# If directly plot each patient directly
plot_clonal_expansion(i=1,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/",self_color=consistent_col) # C1061121
plot_clonal_expansion(i=2,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C107625
plot_clonal_expansion(i=3,lineage_match,metadata,snrna_obj,tumor_anno,percent=99,plot_path="plots/lineages/fishplots/") # C15498
plot_clonal_expansion(i=4,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C176874
plot_clonal_expansion(i=5,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C2542041
plot_clonal_expansion(i=6,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C2751264
plot_clonal_expansion(i=7,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C547104
# only 1 timepoint for C70848
#plot_clonal_expansion(i=8,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C70848
plot_clonal_expansion(i=9,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") # C714384
plot_clonal_expansion(i=10,lineage_match,metadata,snrna_obj,tumor_anno,percent=98,plot_path="plots/lineages/fishplots/") #  C799746
plot_clonal_expansion(i=11,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") #  C1027419
plot_clonal_expansion(i=12,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") #  C1037505
plot_clonal_expansion(i=13,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") #  C2399853
plot_clonal_expansion(i=14,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") #  C2542533
#plot_clonal_expansion(i=15,lineage_match,metadata,snrna_obj,tumor_anno,plot_path="plots/lineages/fishplots/") #  C34809

# plot all patients results into one pdf as extended file
g_list=list()
frac_list = list()
percent=99.5 # default setting
for(i in 1:dim(lineage_match)[1]){
  if(i==3){
    percent=99
  }else if(i==8){
    percent=98
  }
  obj = plot_clonal_expansion(i=i,lineage_match,metadata,
                            snrna_obj,tumor_anno,percent=percent,
                            plot=F,return=T,vlabSize=3,self_color =consistent_col)
  g_list[[i]] <- obj[[1]]  # plot
  frac_list[[i]] <- obj[[2]] # fraction of clones over time
}

pdf(paste0("plots/lineages/fishplots/all_patients_dolphinplot.pdf"),width=25,height=18)
#grid.arrange(grobs=g_list, nrow=7)
grid.arrange(grobs= lapply(g_list, "+", theme(plot.margin=margin(30,0,30,0))),nrow=5)
dev.off()

# By setting certain automated threshold
clones=c()
for(i in 1:length(lineage_match$Patient.ID)){
  fractable=frac_list[[i]]
  # 1. requiring increasing over time
  increase_idx = fractable[,dim(fractable)[2]] >fractable[,1] 
  # 2. require final fraction passing certain threshold
  thresh = 10
  pass_thresh_idx = fractable[,dim(fractable)[2]] > thresh
  remained_clones = rownames(fractable)[increase_idx & pass_thresh_idx]
  clones <- append(clones,paste0(lineage_match$Patient.ID[i],".c",remained_clones))
}
# save expanded clones
saveRDS(clones,"plots/lineages/fishplots/expanded_clones_14patients.rds")

# Or manually Define clonally expanded clones by checking fishplots in each patient
clones = c(paste0(lineage_match$Patient.ID[1],".c",c("2","7")), # C1061121
           paste0(lineage_match$Patient.ID[2],".c",c("1","112")), #C107625
           paste0(lineage_match$Patient.ID[3],".c",c("3","4","170")), # C15498
           paste0(lineage_match$Patient.ID[4],".c",c("1","3")), #C176874
           paste0(lineage_match$Patient.ID[5],".c",c("1","4","200","201","206")), #C2542041
           paste0(lineage_match$Patient.ID[6],".c",c("1","89","114")),#C2751264
           paste0(lineage_match$Patient.ID[7],".c",c("123")), # C547104
           paste0(lineage_match$Patient.ID[8],".c",c("2","5")), # C714384
           paste0(lineage_match$Patient.ID[9],".c",c("2","3")),# C799746
           paste0(lineage_match$Patient.ID[10],".c",c("1","4","126")), # C1027419
           paste0(lineage_match$Patient.ID[11],".c",c("4","9")), # C1037505
           paste0(lineage_match$Patient.ID[12],".c",c("1","9","31","32")), # C2399853
           paste0(lineage_match$Patient.ID[13],".c",c("3","6")), # C2542533
           paste0(lineage_match$Patient.ID[14],".c",c("1","114"))) # C34809
saveRDS(clones,"plots/lineages/fishplots/expanded_clones_14patients_manual_defined.rds")

