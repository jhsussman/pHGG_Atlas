library(Clonalscope)
library(Seurat)
library(ggplot2)
library(stringr)
library(gridExtra)
library(grid)
library("cowplot")
library(dplyr)
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

# load clonally expanded clones from fishplots and tumor clone fractions
clones = readRDS("draft_results/lineages/fishplots/expanded_clones_14patients.rds")

# (1) Integratively run DEG for each patient's tumor cells' all clones
# DEG of tumor vs tumor cells
lineage_match=read.csv("data/new/sample_match/lineage_match_progressiveSeed_14patient.txt",sep="\t",header=T)
colnames(lineage_match) = c("Patient.ID","snRNA_seed","snRNA_fileAccessions","Bulk_sample")
lineage_match = lineage_match[c(1:7,9:15),]
for (i in c(1:dim(lineage_match)[1])){
  patient_id = lineage_match$Patient.ID[i]
  print(paste0(i,", patientID: ",patient_id ))
  seed_snRNA = lineage_match$snRNA_seed[i]
  snRNA_lineages = lineage_match$snRNA_fileAccessions[i]
  bulk = lineage_match$Bulk_sample[i]
  snRNAs = c(seed_snRNA,unlist(str_split(snRNA_lineages,",")))
  timepoints = metadata[snRNAs,c("timepoint")]
  regions = metadata[snRNAs,c("region")]
  
  plot_df_celltype= as.data.frame(cbind(barcodes = names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs]),
                                        snRNA_accession = snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs],
                                        celltype=snrna_obj$cellType[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]))
  plot_df_celltype$Clonalscope_cluster="Filtered/Unknown"
  plot_df_celltype$Neftel = snrna_obj$neftel_type[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]
  plot_df_celltype$TumorAnno = snrna_obj$TumorNormal[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]
  tumor_subclone_bc = intersect(rownames(tumor_anno),
                                names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs]))# then tumor subclones
  plot_df_celltype[tumor_subclone_bc,"TumorAnno"] = tumor_anno[tumor_subclone_bc,"cell_state"]
  
  for (j in 1:length(snRNAs)){
    # load Clonalscope clustering result
    f = snRNAs[j]
    if(j==1){
      #Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly2/",f,"/Cov_obj.rds"))
      Cov0_clusters = levels(as.factor(as.numeric(Cov_obj$result_final$result$Zest)))
    }else{
      #Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly/",patient_id,"/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly_allpatients/",patient_id,"/",f,"/Cov_obj.rds"))
    }
    # load cluster information/annotation
    plot_df_celltype[names(Cov_obj$result_final$result$Zest),"Clonalscope_cluster"] = Cov_obj$result_final$result$Zest
  }
  plot_df_celltype = plot_df_celltype[tumor_subclone_bc,]
  plot_df_celltype$Clonalscope_cluster[is.na(plot_df_celltype$Clonalscope_cluster)] = "Filtered"
  
  # DEG analysis for each patient's tumor cells
  snRNA_sub =  subset(snrna_obj, cells = tumor_subclone_bc)
  new_ident=rep("Unknown/Filtered",length(snRNA_sub@active.ident))
  names(new_ident)= names(snRNA_sub@active.ident)
  new_ident[rownames(plot_df_celltype)] = plot_df_celltype$Clonalscope_cluster
  snRNA_sub@active.ident = as.factor(new_ident) # resetting cluster as Clonalscope's 
  all.markers <- FindAllMarkers(object = snRNA_sub)
  #dir_path = paste0("plots/lineages/interesting_patients/Tumor/",patient_id)
  dir_path = paste0("plots/lineages/DEG/",patient_id)
  dir.create(dir_path,recursive = T)
  #write.table(all.markers,paste0("plots/lineages/interesting_patients/Tumor/",patient_id,"/","integrated_Clonal_DEG.txt"),quote=F,sep="\t")
  write.table(all.markers,paste0("plots/lineages/DEG/",patient_id,"/","integrated_Clonal_DEG.txt"),quote=F,sep="\t")
}

# Tumor vs all cells
for (i in c(1:dim(lineage_match)[1])){
  patient_id = lineage_match$Patient.ID[i]
  print(paste0(i,", patientID: ",patient_id ))
  seed_snRNA = lineage_match$snRNA_seed[i]
  snRNA_lineages = lineage_match$snRNA_fileAccessions[i]
  bulk = lineage_match$Bulk_sample[i]
  snRNAs = c(seed_snRNA,unlist(str_split(snRNA_lineages,",")))
  timepoints = metadata[snRNAs,c("timepoint")]
  regions = metadata[snRNAs,c("region")]
  
  plot_df_celltype= as.data.frame(cbind(barcodes = names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs]),
                                        snRNA_accession = snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs],
                                        celltype=snrna_obj$cellType[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]))
  plot_df_celltype$Clonalscope_cluster="Filtered/Unknown"
  plot_df_celltype$Neftel = snrna_obj$neftel_type[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]
  plot_df_celltype$TumorAnno = snrna_obj$TumorNormal[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]
  tumor_subclone_bc = intersect(rownames(tumor_anno),
                                names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs]))# then tumor subclones
  plot_df_celltype[tumor_subclone_bc,"TumorAnno"] = tumor_anno[tumor_subclone_bc,"cell_state"]
  
  for (j in 1:length(snRNAs)){
    # load Clonalscope clustering result
    f = snRNAs[j]
    if(j==1){
      #Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly2/",f,"/Cov_obj.rds"))
      Cov0_clusters = levels(as.factor(as.numeric(Cov_obj$result_final$result$Zest)))
    }else{
      #Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly/",patient_id,"/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly_allpatients/",patient_id,"/",f,"/Cov_obj.rds"))
    }
    # load cluster information/annotation
    plot_df_celltype[names(Cov_obj$result_final$result$Zest),"Clonalscope_cluster"] = Cov_obj$result_final$result$Zest
  }
  plot_df_celltype = plot_df_celltype[tumor_subclone_bc,]
  plot_df_celltype$Clonalscope_cluster[is.na(plot_df_celltype$Clonalscope_cluster)] = "Filtered"
  
  # DEG analysis for each patient's tumor cells
  snRNA_sub =  subset(snrna_obj, cells = tumor_subclone_bc)
  new_ident=rep("Unknown/Filtered",length(snRNA_sub@active.ident))
  names(new_ident)= names(snRNA_sub@active.ident)
  new_ident[rownames(plot_df_celltype)] = plot_df_celltype$Clonalscope_cluster
  snRNA_sub@active.ident = as.factor(new_ident) # resetting cluster as Clonalscope's 
  all.markers <- FindAllMarkers(object = snRNA_sub)
  #dir_path = paste0("plots/lineages/interesting_patients/Tumor/",patient_id)
  dir_path = paste0("plots/lineages/DEG/",patient_id)
  dir.create(dir_path,recursive = T)
  #write.table(all.markers,paste0("plots/lineages/interesting_patients/Tumor/",patient_id,"/","integrated_Clonal_DEG.txt"),quote=F,sep="\t")
  write.table(all.markers,paste0("plots/lineages/DEG/",patient_id,"/","integrated_Clonal_DEG.txt"),quote=F,sep="\t")
}


# Tumor vs only normal cells




# (2) Plot marker gene/ DEGs of the expand clones in each patient (total 14 working patients)
#reload all DEG results of each patient
deg_mat = c()
for (i in c(1:dim(lineage_match)[1])){
  patient_id = lineage_match$Patient.ID[i]
  print(paste0(i,", patientID: ",patient_id ))
  patient_degs = read.table(paste0("plots/lineages/all_patients/",patient_id,"/","integrated_Clonal_DEG.txt"),sep="\t")
  patient_degs$clone_name = paste0(patient_id,".c",patient_degs$cluster)
  deg_mat <- rbind(deg_mat,patient_degs)
  
  #DEGs[[patient_id]] <- patient_degs
  #sig_clonal_gene_list[[patient_id]] <- patient_degs[abs(patient_degs$avg_log2FC) > 0.5,]
}
deg_mat$size= abs(deg_mat$avg_log2FC)

# 1. plot the most frequenctly occuring DEGs
# filter markers based on criteria
c=as.data.frame(deg_mat) %>% dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(clone_name %in% clones)
# most frequent occuring 20 markers
markers = names(sort(table(c$gene),decreasing = T)[1:20])
c = c%>% dplyr::filter(gene %in% markers)
for (marker in markers){
  for(clone in clones){
    if(!(any(c$clone_name == clone) & any(c$gene == marker))){
      new_row =c(0,0,0,0,1,"Added",marker,clone,0)
      names(new_row) = colnames(c)
      c <- rbind(c,new_row)
    }
  }
}
c$avg_log2FC=as.numeric(c$avg_log2FC);c$size=as.numeric(c$size)
c$gene <- factor(c$gene,levels=markers[length(markers):1]) # ordering levels for plotting
# Plot top 20 occured DEGs
pdf("plots/lineages/DEG/top20_occurance_DEGs.pdf",height=10,width=15)
c%>%ggplot(aes(x=clone_name, y = gene, color = avg_log2FC, size = size)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(
    colours = c("blue","white","red"), #viridis::viridis(20), 
    limits = c(-0.5,0.5), 
    oob = scales::squish, name = 'log2 fold change')
dev.off()

# 2. Plot the top marker of each clone
c=as.data.frame(deg_mat) %>% dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(clone_name %in% clones)
markers = unique((c %>% group_by(clone_name) %>% top_n(n=1,avg_log2FC))$gene)
#markers = names(sort(table(c$gene),decreasing = T)[1:20])
c = as.data.frame(deg_mat) %>% 
  dplyr::filter(clone_name %in% clones) %>% dplyr::filter(gene %in% markers)
for (marker in markers){
  for(clone in clones){
    if(!(any(c$clone_name == clone) & any(c$gene == marker))){
      new_row =c(0,0,0,0,1,"Added",marker,clone,0)
      names(new_row) = colnames(c)
      c <- rbind(c,new_row)
    }
  }
}
c$avg_log2FC=as.numeric(c$avg_log2FC);c$size=as.numeric(c$size)
c$gene <- factor(c$gene,levels=markers[length(markers):1]) # ordering levels for plotting
c$clone_name <- factor(c$clone_name,levels=clones) # ordering levels for plotting
# Plot the literature related DEGs
pdf("plots/lineages/DEG/top_marker_DEG_each_clone.pdf",height=7.5,width=10)
c%>%ggplot(aes(x=clone_name, y = gene, color = avg_log2FC, size = size)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=5)) +
  scale_color_gradientn(
    colours = c("blue","white","red"), #viridis::viridis(20), 
    limits = c(-1,1.5), 
    oob = scales::squish, name = 'log2 fold change') +
  ggtitle("Top Marker Gene of Each Expanded Clone") 
dev.off()

# 3. plot the known CNV related genes OF gbm
markers = c("MDM4","MYCN","PIK3CA","PDGFRA","EGFR","MET","MYC","CCND2","KRAS","CDK4","MDM2", # gains
            "CHD5","CDKN2A/B","PTEN","RB1") # deletions
#markers = names(sort(table(c$gene),decreasing = T)[1:20])
c = as.data.frame(deg_mat) %>% 
  dplyr::filter(clone_name %in% clones) %>% dplyr::filter(gene %in% markers)
for (marker in markers){
  for(clone in clones){
    if(!(any(c$clone_name == clone) & any(c$gene == marker))){
      new_row =c(0,0,0,0,1,"Added",marker,clone,0)
      names(new_row) = colnames(c)
      c <- rbind(c,new_row)
    }
  }
}
c$avg_log2FC=as.numeric(c$avg_log2FC);c$size=as.numeric(c$size)
c$gene <- factor(c$gene,levels=markers[length(markers):1]) # ordering levels for plotting
c$clone_name <- factor(c$clone_name,levels=clones) # ordering levels for plotting
# Plot the literature related DEGs
pdf("plots/lineages/DEG/known_GBM_marker.pdf",height=7.5,width=10)
c%>%ggplot(aes(x=clone_name, y = gene, color = avg_log2FC, size = size)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  theme(axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=5)) +
  scale_color_gradientn(
    colours = c("blue","white","red"), #viridis::viridis(20), 
    #limits = c(-0.5,0.5), 
    oob = scales::squish, name = 'log2 fold change') +
  ggtitle("Top Marker Gene of Each Expanded Clone")
dev.off()

