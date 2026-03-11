library(Clonalscope)
library(Seurat)
library(ggplot2)
library(stringr)
library(gridExtra)
library(grid)
library("cowplot")
library(dplyr)
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)
library(Hmisc)
library(pheatmap)
library(patchwork)
library(ggplotify)
library(grid)

setwd("~/nzhanglab/project/jrong/CPTCA_pHGG/")

source("~/nzhanglab/project/jrong/clonal_scope/scRNACNV/R/Gene_bin_cell_rna_filtered.R")
source("~/nzhanglab/project/jrong/clonal_scope/scRNACNV/R/CreateSegtableNoWGS.R")
source("~/nzhanglab/project/jrong/clonal_scope/scRNACNV/R/Segmentation_bulk.R")

# load data & metadata
snrna_obj <- readRDS("~/nzhanglab/data/CPTCA_Derek_pHGG/Derek_20220731/snRNA_merged_SeuratObj_QCFiltered_DoubletFiltered_downsampled_rPCA_integrated_clustered_withInferCNVTumorNormal_withNeftelTypes.rds")
metadata = read.csv("~/nzhanglab/data/CPTCA_Derek_pHGG/Derek_20220731/HTAN_pHGG_snRNAseq_file_annot_with_Diaz_and_Jabado_samples_2022-07-29.csv")
rownames(metadata) = metadata$fileAccession
#tumor_anno=readRDS("~/nzhanglab/data/CPTCA_Derek_pHGG/Derek_20220731/tumor_metadata_Wenbao.rds")

# read in hg38 genes
gtf=readRDS("~/nzhanglab/project/jrong/clonal_scope/ref/genes_filtered.rds")# GRhC38
gene_id=sapply(strsplit(sapply(strsplit(gtf[,9],";"),"[",1), ' '),'[',2)
gtf$gene_id=gene_id
# read in hg38 chromsomal arms
chrarm_table=read.table("~/nzhanglab/project/jrong/clonal_scope/Clonalscope/data-raw/cytoarm_table_hg38.txt",header=T)
chrarm_table=chrarm_table[order(chrarm_table$chr,chrarm_table$arm,decreasing=F),]
chrarm_table$chr = paste0("chr",chrarm_table$chr)
# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS("~/nzhanglab/project/jrong/clonal_scope/Clonalscope/data-raw/cyclegenes.rds")
# load 200k bp bins for each cell
bin_bed=read.table("~/nzhanglab/project/jrong/clonal_scope/scRNACNV/data-raw/hg38_200kb.windows.bed")
# only keep autosomes (chr 1- 22)
bin_bed=bin_bed[bin_bed[,1] %in% paste0("chr",c(1:22)),]

# for each CNV clone, convert to chromosomal arm level
lineage_match=read.csv("~/nzhanglab/data/CPTCA_Derek_pHGG/Derek_snRNA_20221219/lineage_match_progressiveSeed_14patient.txt",sep="\t",header=T)
colnames(lineage_match) = c("Patient.ID","snRNA_seed","snRNA_fileAccessions","Bulk_sample")
lineage_match = lineage_match[c(1:7,9:15),]

##read clonalscope result
read_scRNACNV_res <- function (scrna,scrna_groups,method="mean"){
  # compute clone average expression
  scrna_clone=NULL
  un=sort(unique(scrna_groups))
  for (clone in un){
    #clone_i <- rowMeans(scrna[,colnames(scrna) %in% gsub("-",".",
    #                                                     names(scrna_groups[scrna_groups == clone]))])
    if(method=='mean'){
      clone_i <- rowMeans(scrna[,colnames(scrna) %in% names(scrna_groups[scrna_groups == clone])],
                          na.rm=T)
    } else if(method=='median'){
      clone_i <- matrixStats::rowMedians(as.matrix(scrna[,colnames(scrna) %in% names(scrna_groups[scrna_groups == clone])]))
    }
    
    scrna_clone <- cbind(scrna_clone, clone_i)
  }
  colnames(scrna_clone) <- un
  
  return(list(scrna_clone=scrna_clone,scrna_groups=scrna_groups,scrna=scrna))
}

## segments/bin to gene
bin_to_arm_gene <-function(cnv_clone_level,seg_table,sub_chroms=NULL,contig="hg38",gtf=NULL, chrarm_table=NULL,method='mean', mode='chrarm'){
  cnv_clone_level=t(cnv_clone_level)
  cnv_gene_level=NULL
  
  
  #if(mode=='infercnv'){
  if(mode == "gene"){
    query_obj = gtf
    query=GenomicRanges::GRanges(gtf[,1], IRanges::IRanges(as.numeric(gtf[,4]),as.numeric(gtf[,5]))) ## cytoband 1-based start and 1-based end
    subject=GenomicRanges::GRanges(paste0('chr',seg_table$chr),IRanges::IRanges(as.numeric(seg_table$start)+1, as.numeric(seg_table$end)))
    ov=findOverlaps(query, subject)
    ov=as.matrix(ov)
    
    for (ii in 1:nrow(query_obj)){
      ## use rtracklayer
      ov_ind = as.numeric(ov[which(ov[,1]==ii),2])
      #chr = paste0("chr",seg_table$chr[idx]); start = seg_table$start[idx]; end=max(seg_table[which(seg_table$chrr==chrr),3, drop=F]);
      if(length(ov_ind)!=0){
        if (method== 'median'){
          a <- apply(cnv_clone_level[ov_ind,, drop=F],2,median)
          
        }else{
          a <- colMeans(cnv_clone_level[ov_ind,, drop=F])
        }
      }else{
        a=rep(NA, ncol(cnv_clone_level))
      }
      cnv_gene_level<- rbind(cnv_gene_level,a)
    }
    rownames(cnv_gene_level) <- gtf$gene_id
    
  }else if(mode == "chrarm"){
    query_obj = chrarm_table
    query=GenomicRanges::GRanges(chrarm_table[,1], IRanges::IRanges(as.numeric(chrarm_table[,3]),as.numeric(chrarm_table[,4]))) ## cytoband 1-based start and 1-based end
    subject=GenomicRanges::GRanges(paste0('chr',seg_table$chr),IRanges::IRanges(as.numeric(seg_table$start)+1, as.numeric(seg_table$end)))
    ov=findOverlaps(query, subject)
    ov=as.matrix(ov)
    
    # for each chromosomal arm
    for (ii in 1:nrow(query_obj)){
      ## use rtracklayer
      ov_ind = as.numeric(ov[which(ov[,1]==ii),2])
      #chr = paste0("chr",seg_table$chr[idx]); start = seg_table$start[idx]; end=max(seg_table[which(seg_table$chrr==chrr),3, drop=F]);
      if(length(ov_ind)!=0){
        if (method== 'median'){
          a <- apply(cnv_clone_level[ov_ind,, drop=F],2,median)
          
        }else if(method== 'mean'){
          a <- colMeans(cnv_clone_level[ov_ind,, drop=F])
        }
      }else{
        a=rep(NA, ncol(cnv_clone_level))
      }
      cnv_gene_level<- rbind(cnv_gene_level,a)
    }
    rownames(cnv_gene_level) <- paste0(chrarm_table[,1],chrarm_table[,2])
  }
  
  return(cnv_gene_level)
}

# convert estimated values to binary gain/loss at arm level
bin_to_arm <-function(cnv_clone_level,seg_table,chrarm_table=NULL,
                      method='binary', mode='chrarm',gain_threshold=1.25,loss_threshold=0.75){
  
  cnv_arm_level=NULL
  gain_mat = NULL
  loss_mat = NULL
  
  
  # find intersecting segments with chromosomal arm
  query_obj = chrarm_table
  query=GenomicRanges::GRanges(chrarm_table[,1], IRanges::IRanges(as.numeric(chrarm_table[,3]),as.numeric(chrarm_table[,4]))) ## cytoband 1-based start and 1-based end
  subject=GenomicRanges::GRanges(paste0('chr',seg_table$chr),IRanges::IRanges(as.numeric(seg_table$start)+1, as.numeric(seg_table$end)))
  ov=findOverlaps(query, subject)
  ov=as.matrix(ov)
  
  # for each chromosomal arm
  for (ii in 1:nrow(query_obj)){
    ## use rtracklayer
    ov_ind = as.numeric(ov[which(ov[,1]==ii),2])
    #chr = paste0("chr",seg_table$chr[idx]); start = seg_table$start[idx]; end=max(seg_table[which(seg_table$chrr==chrr),3, drop=F]);
    if(length(ov_ind)!=0){
      if (method== 'binary'){
        a_gain <-sapply(colnames(cnv_clone_level),function(cluster){
          if(any(cnv_clone_level[ov_ind,as.character(cluster)] > gain_threshold)){
            "gain"
          }else{
            "neu"
          }
        })
        gain_mat = rbind(gain_mat,a_gain)
        # a_gain_length <-sapply(colnames(cnv_clone_level),function(cluster){
        #   seg_idx = which(cnv_clone_level[ov_ind,as.character(cluster)] > gain_threshold)
        #   if(seg_idx >0){
        #     
        #   }else{
        #     0
        #   }
        # })
        a_loss <- sapply(colnames(cnv_clone_level),function(cluster){
          if(any(cnv_clone_level[ov_ind,as.character(cluster)] < loss_threshold)){
            "loss"
          }else{
            "neu"
          }
        })
        loss_mat = rbind(loss_mat,a_loss)
        # a_loss_length <-sapply(colnames(cnv_clone_level),function(cluster){
        #   seg_table_idx = which(cnv_clone_level[ov_ind,as.character(cluster)] < loss_threshold)
        #   if(seg_idx >0){
        # 
        #   }else{
        #     0
        #   }
        # })
      }else if(method== 'mean'){
        a <- colMeans(cnv_clone_level[ov_ind,, drop=F])
      }
    }else{
      a=rep("neu", ncol(cnv_clone_level))
      gain_mat = rbind(gain_mat,a)
      loss_mat = rbind(loss_mat,a)
    }
  }
  rownames(gain_mat) <- paste0(chrarm_table[,1],chrarm_table[,2])
  rownames(loss_mat) <- paste0(chrarm_table[,1],chrarm_table[,2])
  return(list(gain_mat=gain_mat,loss_mat=loss_mat))
}

# initialize patient.clone x chrarm binary mat
chrarm_avg_mat = NULL
# initialize patient.clone x gene mat
chrarm_gene_mat = NULL 
chrarm_gain_mat  = NULL
chrarm_loss_mat = NULL
  
# for each patient
for (i in c(1:dim(lineage_match)[1])){
  patient_id = lineage_match$Patient.ID[i]
  print(paste0(i,", patientID: ",patient_id ))
  seed_snRNA = lineage_match$snRNA_seed[i]
  snRNA_lineages = lineage_match$snRNA_fileAccessions[i]
  bulk = lineage_match$Bulk_sample[i]
  snRNAs = c(seed_snRNA,unlist(str_split(snRNA_lineages,",")))
  #timepoints = metadata[snRNAs,c("timepoint")]
  #regions = metadata[snRNAs,c("region")]
  
  for (j in 1:length(snRNAs)){
    # load Clonalscope clustering result
    f = snRNAs[j]
    if(j==1){
      #Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly2/",f,"/Cov_obj.rds"))
      cnv_profile=Cov_obj$result_final$df_obj$df
      cnv_groups=Cov_obj$result_final$result$Zest
      #Cov0_clusters = levels(as.factor(as.numeric(Cov_obj$result_final$result$Zest)))
    }else{
      #Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly/",patient_id,"/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly_allpatients/",patient_id,"/",f,"/Cov_obj.rds"))
      shared_segments = intersect(colnames(cnv_profile),colnames(Cov_obj$result_final$df_obj$df))
      cnv_profile=rbind(cnv_profile[,shared_segments],Cov_obj$result_final$df_obj$df[,shared_segments])
      cnv_groups= append(cnv_groups,Cov_obj$result_final$result$Zest)
    }
  }
  
  # calculate clone-level CNV average for the patient
  cnv_clone_level = (read_scRNACNV_res(t(cnv_profile),cnv_groups))$scrna_clone
  
  # convert to clone - gene level amplifications (and for correlation calculation)
  seg_table_temp = data.frame(chr=sapply(strsplit(sapply(strsplit(rownames(cnv_clone_level),"-"),"[[",1),"chr"),"[[",2),
                              start=sapply(strsplit(rownames(cnv_clone_level),"-"),"[[",2),
                              end=sapply(strsplit(rownames(cnv_clone_level),"-"),"[[",3))
  
  cnv_clone_gene_level=bin_to_arm_gene(t(cnv_clone_level),
                                       seg_table=seg_table_temp,#Cov_obj$result_final$df_obj$seg_table_filtered,
                                       gtf=gtf,mode="gene")
  colnames(cnv_clone_gene_level) = paste0(patient_id,".c",colnames(cnv_clone_gene_level))
  chrarm_gene_mat = cbind(chrarm_gene_mat,cnv_clone_gene_level)
  
  # convert to clone - average chrarm level cnvs
  cnv_clone_chrarm_level =bin_to_arm_gene(t(cnv_clone_level),
                                          seg_table=seg_table_temp,#Cov_obj$result_final$df_obj$seg_table_filtered,
                                          chrarm_table=chrarm_table)
  colnames(cnv_clone_chrarm_level) = paste0(patient_id,".c",colnames(cnv_clone_chrarm_level))
  chrarm_avg_mat = cbind(chrarm_avg_mat,cnv_clone_chrarm_level)
  
  # convert to clone - binary chrarm level amplifications with set thresholds
  cnv_clone_chrarm_level_binary = bin_to_arm(cnv_clone_level,
                                                  seg_table=seg_table_temp,#Cov_obj$result_final$df_obj$seg_table_filtered,
                                                  chrarm_table=chrarm_table,mode="binary")
  cnv_clone_chrarm_level_gain = cnv_clone_chrarm_level_binary$gain_mat;colnames(cnv_clone_chrarm_level_gain)=paste0(patient_id,".c",colnames(cnv_clone_chrarm_level_gain))
  cnv_clone_chrarm_level_loss = cnv_clone_chrarm_level_binary$loss_mat;colnames(cnv_clone_chrarm_level_loss)=paste0(patient_id,".c",colnames(cnv_clone_chrarm_level_loss))
  
  chrarm_gain_mat = cbind(chrarm_gain_mat, cnv_clone_chrarm_level_gain)
  chrarm_loss_mat = cbind(chrarm_loss_mat, cnv_clone_chrarm_level_loss)
  # call Clonalscope HMM segmentation for each clone - get binary loss/gain and also segment length
  # set new analysis folder
  #save_path = paste0("cnv_clone_level/",patient_id,"/");dir.create(save_path)
  #dir_path <- save_path

  # # reload patient matrix
  # mtx = snrna_obj@assays$RNA@counts[,snrna_obj$patient_id == patient_id]
  # barcodes = as.data.frame(colnames(mtx))
  # features = as.data.frame(unlist(lapply(rownames(mtx),  sub, pattern = "\\.\\d+$", replacement = "")))
  # # convert gene symbol to ensemble id
  # mapIds <- mapIds(EnsDb.Hsapiens.v86, keys = features[,1], keytype = "SYMBOL", columns = c("GENEID"))
  # features[,2] = features[,1] # second column as gene symbol
  # features[,1] = mapIds # first column as ensemble id
  # # filter count matrix
  # rownames(mtx) <- features[,1]
  # mtx <- mtx[!is.na(features[,1]),]
  # features <- features[!is.na(features[,1]),]
  # # remove cycle genes
  # Input_filtered <- FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)
  # # Remove raw inputs
  # rm(mtx); rm(barcodes); rm(features)
  # 
  # bin_mtx=Gen_bin_cell_rna_filtered(bin_bed=bin_bed, barcodes=Input_filtered$barcodes, gene_matrix=Input_filtered$mtx,
  #                                   genes=Input_filtered$features,gtf=gtf)
  # bin_mtx=as.matrix(bin_mtx)
  # saveRDS(bin_mtx,paste0(dir_path,"/","bin_mtx.rds"))
  
  
  #saveRDS(cnv_clone_level,paste0(save_path,"clonalscope_cnv_clone_segment_average.rds"))
  #saveRDS(cnv_clone_chrarm_level,paste0(save_path,"clonalscope_cnv_clone_chrarm_average.rds"))
}

saveRDS(chrarm_gene_mat,"cnv_clone_level/clone_gene_mat.rds")
saveRDS(chrarm_gain_mat,"cnv_clone_level/chrarm_gain_mat.rds")
saveRDS(chrarm_loss_mat,"cnv_clone_level/chrarm_loss_mat.rds")
saveRDS(chrarm_avg_mat,"cnv_clone_level/chrarm_avg_mat.rds")

###### Reload results for plotting 
chrarm_gene_mat=readRDS("cnv_clone_level/clone_gene_mat.rds")
chrarm_gain_mat=readRDS("cnv_clone_level/chrarm_gain_mat.rds")
chrarm_loss_mat=readRDS("cnv_clone_level/chrarm_loss_mat.rds")
chrarm_avg_mat=readRDS("cnv_clone_level/chrarm_avg_mat.rds")

#### plot overall summary
# correlation at gene level
gene_cor = rcorr(chrarm_gene_mat)
# correlation at chromosomal arm level
chrarm_avg_cor = rcorr(chrarm_avg_mat)
# correlation at binary chromosomal arm level
binary_gain_loss_mat = rbind(chrarm_gain_mat,chrarm_loss_mat)
binary_gain_loss_mat[binary_gain_loss_mat == "neu"] = 0
binary_gain_loss_mat[binary_gain_loss_mat == "loss"] = -1
binary_gain_loss_mat[binary_gain_loss_mat == "gain"] = 1
chrarm_binary_cor = rcorr(binary_gain_loss_mat,type="pearson")
chrarm_binary_cor$r[is.na(chrarm_binary_cor$r)] = 0

col <- colorRampPalette(c("blue","white","red"))(100)
color_breaks = c(seq(-1,1,length.out=101))

# add in annotation of clones
expanded_clones = readRDS("draft_results/lineages/fishplots/expanded_clones_14patients.rds")
clone_annotation = data.frame(clones = colnames(chrarm_gene_mat))
clone_annotation$patient = factor(sapply(strsplit(clone_annotation$clones,".c"),"[[",1))
clone_annotation$expand = "Not expanded";
clone_annotation$expand[clone_annotation$clones %in% expanded_clones] = "Expanded"
clone_annotation$expand = factor(clone_annotation$expand)
clone_colors = c("#27AAe1","#F9ED32","#2CA02C","#DA1C5C","#F26522","#2E3192","#E377C2","#7F7F7F","#006A02")
names(clone_colors) = c("#27AAe1","#F9ED32","#2CA02C","#DA1C5C","#F26522","#2E3192","#E377C2","#7F7F7F","#006A02")
clone_annotation$clone_color = factor(unlist(lapply(unique(clone_annotation$patient),function(pid){clone_colors[1:sum(clone_annotation$patient == pid)]})))
rownames(clone_annotation) = clone_annotation$clones; clone_annotation = clone_annotation[,-1]

# Specify colors
patient_col = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(clone_annotation$patient)))
names(patient_col) = levels(clone_annotation$patient)
expand_cols = c("red", "grey");names(expand_cols)=c(levels(clone_annotation$expand))
ann_colors = list(patient = patient_col, expand = expand_cols,clone_color=clone_colors)

pdf("cnv_clone_level/Correlation_between_all_clones.pdf",height=6,width=6)
pheatmap(gene_cor$r,col=col,breaks = color_breaks,border_color = NA, 
         fontsize=5,symm=T,scale="none",
         annotation_col=clone_annotation,annotation_colors = ann_colors,
         cellwidth=3,cellheight=3,
         main="Pearson Correlation of Clone CNVs at Gene Level")
pheatmap(chrarm_avg_cor$r,col=col,breaks = color_breaks,border_color = NA, fontsize=5,symm=T,scale="none",
         annotation_col=clone_annotation,annotation_colors = ann_colors,
         cellwidth=3,cellheight=3,
         main="Pearson Correlation of Clone Average CNVs at Chromsomal Arm Level")
pheatmap(chrarm_binary_cor$r,col=col,breaks = color_breaks,border_color = NA,fontsize=5,symm=T,scale="none",
         annotation_col=clone_annotation,annotation_colors = ann_colors,
         cellwidth=3,cellheight=3,
         main="Pearson Correlation of Clone Binary CNVs at Chromsomal Arm Level")
dev.off()

# load expanded clones
clones = readRDS("draft_results/lineages/fishplots/expanded_clones_14patients.rds")
# correlation at gene level
gene_cor = rcorr(chrarm_gene_mat[,clones])
# correlation at chromosomal arm level
chrarm_avg_cor = rcorr(chrarm_avg_mat[,clones])
# correlation at binary chromosomal arm level
binary_gain_loss_mat = rbind(chrarm_gain_mat[,clones],chrarm_loss_mat[,clones])
binary_gain_loss_mat[binary_gain_loss_mat == "neu"] = 0
binary_gain_loss_mat[binary_gain_loss_mat == "loss"] = -1
binary_gain_loss_mat[binary_gain_loss_mat == "gain"] = 1
chrarm_binary_cor = rcorr(binary_gain_loss_mat,type="pearson")
chrarm_binary_cor$r[is.na(chrarm_binary_cor$r)] = 0

col <- colorRampPalette(c("blue","white","red"))(100)
color_breaks = c(seq(-1,1,length.out=101))

pdf("cnv_clone_level/Correlation_between_expanded_clones.pdf",height=6,width=6)
pheatmap(gene_cor$r,col=col,breaks = color_breaks,border_color = NA, fontsize=5,symm=T,scale="none",
         annotation_col=clone_annotation,annotation_colors = ann_colors,
         cellwidth=8,cellheight=8,
         main="Pearson Correlation of Clone CNVs at Gene Level")
pheatmap(chrarm_avg_cor$r,col=col,breaks = color_breaks,border_color = NA, fontsize=5,symm=T,scale="none",
         annotation_col=clone_annotation,annotation_colors = ann_colors,
         cellwidth=8,cellheight=8,
         main="Pearson Correlation of Clone Average CNVs at Chromsomal Arm Level")
pheatmap(chrarm_binary_cor$r,col=col,breaks = color_breaks,border_color = NA,fontsize=5,symm=T,scale="none",
         annotation_col=clone_annotation,annotation_colors = ann_colors,
         cellwidth=8,cellheight=8,
         main="Pearson Correlation of Clone Binary CNVs at Chromsomal Arm Level")
dev.off()

# Calculate clone size
clone_table = c()
for (i in c(1:dim(lineage_match)[1])){
  patient_id = lineage_match$Patient.ID[i]
  print(paste0(i,", patientID: ",patient_id ))
  seed_snRNA = lineage_match$snRNA_seed[i]
  snRNA_lineages = lineage_match$snRNA_fileAccessions[i]
  bulk = lineage_match$Bulk_sample[i]
  snRNAs = c(seed_snRNA,unlist(str_split(snRNA_lineages,",")))
  #timepoints = metadata[snRNAs,c("timepoint")]
  #regions = metadata[snRNAs,c("region")]
  
  Zest_list <- c()
  for (j in 1:length(snRNAs)){
    # load Clonalscope clustering result
    f = snRNAs[j]
    if(j==1){
      #Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly2/",f,"/Cov_obj.rds"))
      Cov0_clusters = levels(as.factor(as.numeric(Cov_obj$result_final$result$Zest)))
      Zest_list <- Cov_obj$result_final$result$Zest
    }else{
      #Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly/",patient_id,"/",f,"/Cov_obj.rds"))
      Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly_allpatients/",patient_id,"/",f,"/Cov_obj.rds"))
      Zest_list <- append(Zest_list,Cov_obj$result_final$result$Zest)
    }
  }
  # count clone size
  clone_tab_i = as.data.frame(table(Zest_list)); clone_tab_i$Zest_list = paste0(patient_id,".c",clone_tab_i$Zest_list)
  clone_table = rbind(clone_table,clone_tab_i)
}
rownames(clone_table) = clone_table$Zest_list


# plot gain/loss across all patients clones
col <- colorRampPalette(c("blue","white","red"))(100)
color_breaks = c(seq(0.5,1.5,length.out=101))
chrarm_avg_mat[is.na(chrarm_avg_mat)]= 1 # copy neutral
# plot chr arm x clones average profile
avg_heat <- pheatmap(chrarm_avg_mat,col=col,breaks = color_breaks,border_color = NA,
                     annotation_col=clone_annotation,annotation_colors = ann_colors,
                     cluster_rows=F,cluster_cols=F,main="Average CNV across all clones",fontsize=5)
# plot binary frequency
gainloss_freq = as.data.frame(chrarm_gain_mat); gainloss_freq[gainloss_freq=="neu"] = 0;gainloss_freq[gainloss_freq=="gain"] = 1;
gainloss_freq=sapply(gainloss_freq,as.numeric);rownames(gainloss_freq) = rownames(chrarm_gain_mat)
gainloss_freq=data.frame(chrarm=factor(rownames(gainloss_freq),levels=rev(rownames(gainloss_freq))),freq=rowSums(gainloss_freq),gain_loss="Gain")
gainloss_freq$binary_pval = NA;gainloss_freq$binary_padj = NA;
gainloss_freq$cellsize_pval = NA;gainloss_freq$cellsize_padj = NA;
# Test for significance of the gain/loss frequency
total_clone_num = dim(chrarm_gain_mat)[2]
total_cell_num = sum(clone_table$Freq)
for(i in 1:dim(gainloss_freq)[1]){
  # contingency test based on binary frequency
  obs_table = Matrix(c(gainloss_freq[i,"freq"],total_clone_num-gainloss_freq[i,"freq"], # in this segment, gain and not gain freq
                       sum(gainloss_freq$freq)-gainloss_freq[i,"freq"], # rest - gain freq
                       total_clone_num*dim(gainloss_freq)[1]- sum(gainloss_freq$freq)-gainloss_freq[i,"freq"]), # rest- not gain freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("gain","not gain")
  gainloss_freq[i,"binary_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
  # contingency test based on cell-size frequency
  rest_seg_not_gain = sum(sapply(1:dim(gainloss_freq)[1],function(j){sum(clone_table[chrarm_gain_mat[j,] == "gain","Freq"])})) - sum(clone_table[chrarm_gain_mat[i,] == "gain","Freq"])
  obs_table = Matrix(c(sum(clone_table[chrarm_gain_mat[i,] == "gain","Freq"]), # in this segment -gain
                       total_cell_num-sum(clone_table[chrarm_gain_mat[i,] == "gain","Freq"]), # in this segment -not gain
                       rest_seg_not_gain, # rest segments - gain freq
                       total_cell_num*dim(gainloss_freq)[1]- rest_seg_not_gain), # rest- not gain freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("gain","not gain")
  gainloss_freq[i,"cellsize_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
}
gainloss_freq$binary_padj = p.adjust(gainloss_freq$binary_pval,method = "fdr");
gainloss_freq$cellsize_padj = p.adjust(gainloss_freq$cellsize_pval,method = "fdr");

loss_freq=as.data.frame(chrarm_loss_mat); loss_freq[loss_freq=="neu"] = 0; loss_freq[loss_freq=="loss"] = 1;
loss_freq=sapply(loss_freq,as.numeric);rownames(loss_freq) = rownames(chrarm_loss_mat)
loss_freq=data.frame(chrarm=factor(rownames(loss_freq),levels=rev(rownames(loss_freq))),freq=rowSums(loss_freq),gain_loss="Loss")
loss_freq$binary_pval = NA;loss_freq$binary_padj = NA;
loss_freq$cellsize_pval = NA;loss_freq$cellsize_padj = NA;
# Test for significance of the gain/loss frequency
total_clone_num = dim(chrarm_loss_mat)[2]
total_cell_num = sum(clone_table$Freq)
for(i in 1:dim(loss_freq)[1]){
  # contingency test based on binary frequency
  obs_table = Matrix(c(loss_freq[i,"freq"],total_clone_num-loss_freq[i,"freq"], # in this segment, gain and not gain freq
                       sum(loss_freq$freq)-loss_freq[i,"freq"], # rest - gain freq
                       total_clone_num*dim(loss_freq)[1]- sum(loss_freq$freq)-loss_freq[i,"freq"]), # rest- not gain freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("loss","not loss")
  loss_freq[i,"binary_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
  # contingency test based on cell-size frequency
  rest_seg_not_gain = sum(sapply(1:dim(loss_freq)[1],function(j){sum(clone_table[chrarm_loss_mat[j,] == "loss","Freq"])})) - sum(clone_table[chrarm_gain_mat[i,] == "gain","Freq"])
  obs_table = Matrix(c(sum(clone_table[chrarm_loss_mat[i,] == "loss","Freq"]), # in this segment -loss
                       total_cell_num-sum(clone_table[chrarm_loss_mat[i,] == "loss","Freq"]), # in this segment -not loss
                       rest_seg_not_gain, # rest segments - loss freq
                       total_cell_num*dim(gainloss_freq)[1]- rest_seg_not_gain), # rest- not loss freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("loss","not loss")
  loss_freq[i,"cellsize_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
}
loss_freq$binary_padj = p.adjust(loss_freq$binary_pval,method = "fdr");
loss_freq$cellsize_padj = p.adjust(loss_freq$cellsize_pval,method = "fdr");

gainloss_freq = rbind(gainloss_freq,loss_freq)
gainloss_freq$ratio = gainloss_freq$freq/dim(chrarm_gain_mat)[2]
saveRDS(gainloss_freq,"cnv_clone_level/gainloss_freq_all_clones.rds")

gainloss_freq$label = -log10(gainloss_freq$binary_padj)
gainloss_freq$label = unlist(lapply(gainloss_freq$label,function(x){
  if(x >=3){ # p < 1e-3
    "***"
  }else if(x >=2){ # p < 1e-2
    "**"
  }else if(x >=-log10(0.05)){# p < 0.05
    "*"
  }else{
    NA
  }
}))
gainloss_freq$label_position = gainloss_freq$ratio + 0.08

gl_color = c(col[1],col[100]);names(gl_color) = c("Loss","Gain")
bar <- ggplot(gainloss_freq,aes(y=chrarm,x=ratio,fill=gain_loss)) + geom_bar(stat="identity") + facet_wrap(~gain_loss) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     text=element_text(size=5),)+ 
  xlab("Ratio") + ylab("") + xlim(0, max(gainloss_freq$ratio)+0.15)+
  scale_fill_manual(values=gl_color) + labs(fill="Gain/Loss") +
  geom_text(aes(label = label, x=label_position,y=chrarm))
bar
  
# plot chr arm x clones binary gain/loss
binary_mat = chrarm_gain_mat
binary_mat[(chrarm_gain_mat == "neu") & (chrarm_loss_mat == "neu")] = "Neutral";
binary_mat[chrarm_gain_mat == "gain"] = "Gain"
binary_mat[chrarm_loss_mat == "loss"] = "Loss"
binary_mat[(chrarm_gain_mat == "gain") & (chrarm_loss_mat == "loss")] = "Gain & Loss"
rownames(binary_mat) = rownames(chrarm_gain_mat);colnames(binary_mat) = colnames(chrarm_gain_mat);
long_binary_df=data.frame(chrarm=rownames(binary_mat),clones=rep(colnames(binary_mat),each=length(rownames(binary_mat))),
                          value=c(binary_mat))
long_binary_df$chrarm = factor(long_binary_df$chrarm,levels=rev(rownames(binary_mat)))
long_binary_df$clones=factor(long_binary_df$clones,levels=colnames(binary_mat))
#long_binary_df$value = factor(long_binary_df$value,levels=c("Gain","Neutral","Loss","Gain & Loss"))
binary_color = c("red","white","blue","purple");names(binary_color) = c("Gain","Neutral","Loss","Gain & Loss") 
#binary_heat = pheatmap(as.factor(binary_mat),border_color = NA,
#                       cluster_rows=F,cluster_cols=F,main="Binary CNV Gain and Loss across all clones")
binary_heat = ggplot(long_binary_df,aes(x=clones,y=chrarm))+geom_tile(aes(fill=value),color="grey") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        text=element_text(size=5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_fill_manual(values=binary_color) + 
  ggtitle("Binary CNV Gain and Loss aross all clones") + labs(fill="Gain/Loss") + coord_fixed()

pdf("cnv_clone_level/cnv_summary_heatmap_across_all_clones.pdf",width=20,height=12)
plot_grid(as.ggplot(avg_heat),bar, nrow = 1,rel_widths = c(3/4,1/4));
plot_grid(binary_heat,bar, nrow = 1,rel_widths = c(5/8,3/8));
dev.off()

# readjust for plotting
# plot chr arm x clones average profile
# patient break in columns
col_gap = sort(unlist(lapply(unique(clone_annotation$patient),
                function(i){which(clone_annotation$patient == unique(clone_annotation$patient)[i])[length(which(clone_annotation$patient == unique(clone_annotation$patient)[i]))]})))
avg_heat_gapped <- pheatmap(chrarm_avg_mat,col=col,breaks = color_breaks,border_color = NA,
                     annotation_col=clone_annotation,annotation_colors = ann_colors,
                     cellwidth=4,cellheight=2,
                     width=5,height=2.5,
                     gaps_col = col_gap,
                     cluster_rows=F,cluster_cols=F,main="Average CNV across all clones",fontsize=5)
avg_heat <- pheatmap(chrarm_avg_mat,col=col,breaks = color_breaks,border_color = NA,
                     annotation_col=clone_annotation,annotation_colors = ann_colors,
                     cellwidth=4,cellheight=2,
                     width=5,height=2.5,
                     cluster_rows=F,cluster_cols=F,main="Average CNV across all clones",fontsize=5)
bar <- ggplot(gainloss_freq,aes(y=chrarm,x=ratio,fill=gain_loss)) + geom_bar(stat="identity") + facet_wrap(~gain_loss) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     text=element_text(size=5),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+ 
  xlab("Ratio") + ylab("") + xlim(0, max(gainloss_freq$ratio)+0.15)+
  scale_fill_manual(values=gl_color) + labs(fill="Gain/Loss") +
  geom_text(aes(label = label, x=label_position,y=chrarm)) 
bar
pdf("cnv_clone_level/cnv_summary_sized.pdf",width=10,height=3)
plot_grid(as.ggplot(avg_heat_gapped),bar, nrow = 1,rel_widths = c(7/8,1/8));
plot_grid(as.ggplot(avg_heat),bar, nrow = 1,rel_widths = c(7/8,1/8));
plot_grid(binary_heat,bar, nrow = 1,rel_widths = c(5/8,3/8));
dev.off()

## plot gain/loss across in expanded clones
expanded_clones = readRDS("draft_results/lineages/fishplots/expanded_clones_14patients.rds")
col <- colorRampPalette(c("blue","white","red"))(100)
color_breaks = c(seq(0.5,1.5,length.out=101))
chrarm_avg_mat[is.na(chrarm_avg_mat)]= 1 # copy neutral
# plot chr arm x clones average profile
avg_heat <- pheatmap(chrarm_avg_mat[,expanded_clones],col=col,breaks = color_breaks,
                     border_color = NA,cluster_rows=F,cluster_cols=F,main="Average CNV across expanded clones",
                     annotation_col=clone_annotation,annotation_colors = ann_colors,
                     fontsize=10)
# plot binary frequency
gainloss_freq = as.data.frame(chrarm_gain_mat[,expanded_clones]); gainloss_freq[gainloss_freq=="neu"] = 0;gainloss_freq[gainloss_freq=="gain"] = 1;
gainloss_freq=sapply(gainloss_freq,as.numeric);rownames(gainloss_freq) = rownames(chrarm_gain_mat[,expanded_clones])
gainloss_freq=data.frame(chrarm=factor(rownames(gainloss_freq),levels=rev(rownames(gainloss_freq))),freq=rowSums(gainloss_freq),gain_loss="Gain")
gainloss_freq$binary_pval = NA;gainloss_freq$binary_padj = NA;
gainloss_freq$cellsize_pval = NA;gainloss_freq$cellsize_padj = NA;
# Test for significance of the gain/loss frequency
total_clone_num = dim(chrarm_gain_mat[,expanded_clones])[2]
total_cell_num = sum(clone_table$Freq)
for(i in 1:dim(gainloss_freq)[1]){
  # contingency test based on binary frequency
  obs_table = Matrix(c(gainloss_freq[i,"freq"],total_clone_num-gainloss_freq[i,"freq"], # in this segment, gain and not gain freq
                       sum(gainloss_freq$freq)-gainloss_freq[i,"freq"], # rest - gain freq
                       total_clone_num*dim(gainloss_freq)[1]- sum(gainloss_freq$freq)-gainloss_freq[i,"freq"]), # rest- not gain freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("gain","not gain")
  gainloss_freq[i,"binary_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
  # contingency test based on cell-size frequency
  rest_seg_not_gain = sum(sapply(1:dim(gainloss_freq)[1],function(j){sum(clone_table[chrarm_gain_mat[,expanded_clones][j,] == "gain","Freq"])})) - sum(clone_table[chrarm_gain_mat[,expanded_clones][i,] == "gain","Freq"])
  obs_table = Matrix(c(sum(clone_table[chrarm_gain_mat[,expanded_clones][i,] == "gain","Freq"]), # in this segment -gain
                       total_cell_num-sum(clone_table[chrarm_gain_mat[,expanded_clones][i,] == "gain","Freq"]), # in this segment -not gain
                       rest_seg_not_gain, # rest segments - gain freq
                       total_cell_num*dim(gainloss_freq)[1]- rest_seg_not_gain), # rest- not gain freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("gain","not gain")
  gainloss_freq[i,"cellsize_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
}
gainloss_freq$binary_padj = p.adjust(gainloss_freq$binary_pval,method = "fdr");
gainloss_freq$cellsize_padj = p.adjust(gainloss_freq$cellsize_pval,method = "fdr");

loss_freq=as.data.frame(chrarm_loss_mat[,expanded_clones]); loss_freq[loss_freq=="neu"] = 0; loss_freq[loss_freq=="loss"] = 1;
loss_freq=sapply(loss_freq,as.numeric);rownames(loss_freq) = rownames(chrarm_loss_mat[,expanded_clones])
loss_freq=data.frame(chrarm=factor(rownames(loss_freq),levels=rev(rownames(loss_freq))),freq=rowSums(loss_freq),gain_loss="Loss")
loss_freq$binary_pval = NA;loss_freq$binary_padj = NA;
loss_freq$cellsize_pval = NA;loss_freq$cellsize_padj = NA;
# Test for significance of the gain/loss frequency
total_clone_num = dim(chrarm_loss_mat[,expanded_clones])[2]
total_cell_num = sum(clone_table$Freq)
for(i in 1:dim(loss_freq)[1]){
  # contingency test based on binary frequency
  obs_table = Matrix(c(loss_freq[i,"freq"],total_clone_num-loss_freq[i,"freq"], # in this segment, gain and not gain freq
                       sum(loss_freq$freq)-loss_freq[i,"freq"], # rest - gain freq
                       total_clone_num*dim(loss_freq)[1]- sum(loss_freq$freq)-loss_freq[i,"freq"]), # rest- not gain freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("loss","not loss")
  loss_freq[i,"binary_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
  # contingency test based on cell-size frequency
  rest_seg_not_gain = sum(sapply(1:dim(loss_freq)[1],function(j){sum(clone_table[chrarm_loss_mat[,expanded_clones][j,] == "loss","Freq"])})) - sum(clone_table[chrarm_gain_mat[,expanded_clones][i,] == "gain","Freq"])
  obs_table = Matrix(c(sum(clone_table[chrarm_loss_mat[i,] == "loss","Freq"]), # in this segment -loss
                       total_cell_num-sum(clone_table[chrarm_loss_mat[,expanded_clones][i,] == "loss","Freq"]), # in this segment -not loss
                       rest_seg_not_gain, # rest segments - loss freq
                       total_cell_num*dim(gainloss_freq)[1]- rest_seg_not_gain), # rest- not loss freq
                     nrow=2,byrow=T)
  rownames(obs_table) =c("this segment","rest segments");colnames(obs_table) = c("loss","not loss")
  loss_freq[i,"cellsize_pval"] = (fisher.test(as.matrix(obs_table),alternative="greater"))$p.value
  
}
loss_freq$binary_padj = p.adjust(loss_freq$binary_pval,method = "fdr");
loss_freq$cellsize_padj = p.adjust(loss_freq$cellsize_pval,method = "fdr");

gainloss_freq = rbind(gainloss_freq,loss_freq)
gainloss_freq$ratio = gainloss_freq$freq/dim(chrarm_gain_mat[,expanded_clones])[2]

saveRDS(gainloss_freq,"cnv_clone_level/gainloss_freq_expanded_clones.rds")

gl_color = c(col[1],col[100]);names(gl_color) = c("Loss","Gain")
gainloss_freq$label = -log10(gainloss_freq$binary_padj)
gainloss_freq$label = unlist(lapply(gainloss_freq$label,function(x){
  if(x >=3){ # p < 1e-3
    "***"
  }else if(x >=2){ # p < 1e-2
    "**"
  }else if(x >=-log10(0.05)){# p < 0.05
    "*"
  }else{
    NA
  }
}))
gainloss_freq$label_position = gainloss_freq$ratio + 0.08
bar <- ggplot(gainloss_freq,aes(y=chrarm,x=ratio,fill=gain_loss)) + geom_bar(stat="identity") + facet_wrap(~gain_loss) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     text=element_text(size=10),)+ 
  xlab("Ratio") + ylab("") + xlim(0, max(gainloss_freq$ratio)+0.15)+
  scale_fill_manual(values=gl_color) + labs(fill="Gain/Loss") +
  geom_text(aes(label = label, x=label_position,y=chrarm))
bar
# plot chr arm x clones binary gain/loss
binary_mat = chrarm_gain_mat[,expanded_clones]
binary_mat[(chrarm_gain_mat[,expanded_clones] == "neu") & (chrarm_loss_mat[,expanded_clones] == "neu")] = "Neutral";
binary_mat[chrarm_gain_mat[,expanded_clones] == "gain"] = "Gain"
binary_mat[chrarm_loss_mat[,expanded_clones] == "loss"] = "Loss"
binary_mat[(chrarm_gain_mat[,expanded_clones] == "gain") & (chrarm_loss_mat[,expanded_clones] == "loss")] = "Gain & Loss"
rownames(binary_mat) = rownames(chrarm_gain_mat[,expanded_clones]);colnames(binary_mat) = colnames(chrarm_gain_mat[,expanded_clones]);
long_binary_df=data.frame(chrarm=rownames(binary_mat),clones=rep(colnames(binary_mat),each=length(rownames(binary_mat))),
                          value=c(binary_mat))
long_binary_df$chrarm = factor(long_binary_df$chrarm,levels=rev(rownames(binary_mat)))
long_binary_df$clones=factor(long_binary_df$clones,levels=colnames(binary_mat))
#long_binary_df$value = factor(long_binary_df$value,levels=c("Gain","Neutral","Loss","Gain & Loss"))
binary_color = c("red","white","blue","purple");names(binary_color) = c("Gain","Neutral","Loss","Gain & Loss") 
#binary_heat = pheatmap(as.factor(binary_mat),border_color = NA,
#                       cluster_rows=F,cluster_cols=F,main="Binary CNV Gain and Loss across all clones")
binary_heat = ggplot(long_binary_df,aes(x=clones,y=chrarm))+geom_tile(aes(fill=value),color="grey") +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_fill_manual(values=binary_color) + 
  ggtitle("Binary CNV Gain and Loss aross expanded clones") + labs(fill="Gain/Loss") +coord_fixed()

pdf("cnv_clone_level/cnv_summary_heatmap_across_expanded_clones.pdf",width=12,height=6)
plot_grid(as.ggplot(avg_heat),bar, nrow = 1,rel_widths = c(3/4,1/4));
plot_grid(binary_heat,bar, nrow = 1,rel_widths = c(1/2,1/2));
dev.off()

pdf("cnv_clone_level/cnv_summary_heatmap_across_expanded_clones.pdf",width=12,height=6)
plot_grid(as.ggplot(avg_heat),bar, nrow = 1,rel_widths = c(3/4,1/4));
plot_grid(binary_heat,bar, nrow = 1,rel_widths = c(1/2,1/2));
dev.off()



# Test for expanded clones vs non-expanded clones

# 1. plot ratio plot separately
# plot binary frequency by gain, loss, expanded, non-expanded
rest_clones=setdiff(colnames(chrarm_avg_mat),expanded_clones)
# plot binary frequency
gainloss_freq_expand=readRDS("cnv_clone_level/gainloss_freq_expanded_clones.rds")
gainloss_freq_expand$expand="Expanded"

gainloss_freq = as.data.frame(chrarm_gain_mat[,rest_clones]); gainloss_freq[gainloss_freq=="neu"] = 0;gainloss_freq[gainloss_freq=="gain"] = 1;
gainloss_freq=sapply(gainloss_freq,as.numeric);rownames(gainloss_freq) = rownames(chrarm_gain_mat[,rest_clones])
gainloss_freq=data.frame(chrarm=factor(rownames(gainloss_freq),levels=rev(rownames(gainloss_freq))),freq=rowSums(gainloss_freq),gain_loss="Gain")
loss_freq=as.data.frame(chrarm_loss_mat[,rest_clones]); loss_freq[loss_freq=="neu"] = 0; loss_freq[loss_freq=="loss"] = 1;
loss_freq=sapply(loss_freq,as.numeric);rownames(loss_freq) = rownames(chrarm_loss_mat[,rest_clones])
loss_freq=data.frame(chrarm=factor(rownames(loss_freq),levels=rev(rownames(loss_freq))),freq=rowSums(loss_freq),gain_loss="Loss")
gainloss_freq = rbind(gainloss_freq,loss_freq)
gainloss_freq$ratio = gainloss_freq$freq/dim(chrarm_gain_mat[,rest_clones])[2]
gainloss_freq$expand="Non Expanded"

gainloss_freq_all = rbind(gainloss_freq_expand[,colnames(gainloss_freq)],gainloss_freq)
saveRDS(gainloss_freq_all,"cnv_clone_level/gainloss_freq_expanded_non_expanded_clones.rds")

gainloss_freq_all = readRDS("cnv_clone_level/gainloss_freq_expanded_non_expanded_clones.rds")
gainloss_freq_all$binary_pval = NA;gainloss_freq_all$binary_padj=NA;
gainloss_freq_all$cellsize_pval = NA;gainloss_freq_all$cellsize_padj=NA;
total_cell_num = sum(clone_table$Freq)
total_expand_cell = sum(clone_table[expanded_clones,"Freq"])
total_non_expand_cell = sum(clone_table[rest_clones,"Freq"])

for(i in 1:length(unique(gainloss_freq_all$chrarm))){
  seg = unique(gainloss_freq_all$chrarm[i])
  print(seg)
  gain_idx = (gainloss_freq_all$chrarm == seg & gainloss_freq_all$gain_loss == "Gain")
  obs = as.matrix(Matrix(c(gainloss_freq_all[gain_idx,"freq"],
                 c(length(expanded_clones), length(rest_clones))-c(gainloss_freq_all[gain_idx,"freq"])),
               byrow=F,nrow=2))
  rownames(obs)=c("Expanded","Non Expanded");colnames(obs) = c("Gain","Not Gain")
  gainloss_freq_all[gain_idx,"binary_pval"] = (fisher.test(obs,alternative="greater"))$p.value
    #(prop.test(x = gainloss_freq_all[gain_idx,"freq"],n = c(length(expanded_clones), length(rest_clones)),alternative="greater"))$p.value
  loss_idx=(gainloss_freq_all$chrarm == seg & gainloss_freq_all$gain_loss == "Loss")
  obs = as.matrix(Matrix(c(gainloss_freq_all[loss_idx,"freq"],
                           c(length(expanded_clones), length(rest_clones))-c(gainloss_freq_all[loss_idx,"freq"])),
                         byrow=F,nrow=2))
  gainloss_freq_all[loss_idx,"binary_pval"] = (fisher.test(obs,alternative="greater"))$p.value
  #(prop.test(x = gainloss_freq_all[loss_idx,"freq"],n = c(length(expanded_clones), length(rest_clones)),alternative="greater"))$p.value
  # cell size exact table
  obs = as.matrix(Matrix(c(sum(clone_table[expanded_clones,][chrarm_gain_mat[,expanded_clones][i,] == "gain","Freq"]), # in expanded -gain
                       total_expand_cell-sum(clone_table[expanded_clones,][chrarm_gain_mat[,expanded_clones][i,] == "gain","Freq"]), # expand  -not gain
                       sum(clone_table[rest_clones,][chrarm_gain_mat[,rest_clones][i,] == "gain","Freq"]), # in non expanded -gain
                       total_non_expand_cell- sum(clone_table[rest_clones,][chrarm_gain_mat[,rest_clones][i,] == "gain","Freq"])), # non expand- not gain freq
                     nrow=2,byrow=T))
  rownames(obs)=c("Expanded","Non Expanded");colnames(obs) = c("Gain","Not Gain")
  gainloss_freq_all[gain_idx,"cellsize_pval"] = (fisher.test(obs,alternative="greater"))$p.value
  obs = as.matrix(Matrix(c(sum(clone_table[expanded_clones,][chrarm_loss_mat[,expanded_clones][i,] == "loss","Freq"]), # in expanded - loss
                           total_expand_cell-sum(clone_table[expanded_clones,][chrarm_loss_mat[,expanded_clones][i,] == "loss","Freq"]), # expand  -not loss
                           sum(clone_table[rest_clones,][chrarm_loss_mat[,rest_clones][i,] == "gain","Freq"]), # in non expanded - loss
                           total_non_expand_cell- sum(clone_table[rest_clones,][chrarm_loss_mat[,rest_clones][i,] == "loss","Freq"])), # non expand- not loss freq
                         nrow=2,byrow=T))
  rownames(obs)=c("Expanded","Non Expanded");colnames(obs) = c("Loss","Not Loss")
  gainloss_freq_all[loss_idx,"cellsize_pval"] = (fisher.test(obs,alternative="greater"))$p.value
}
gainloss_freq_all[gainloss_freq_all$expand == "Non Expanded","binary_pval"]=NA
gainloss_freq_all$binary_padj[1:88] = p.adjust(gainloss_freq_all$binary_pval[1:88],method="fdr") 
gainloss_freq_all[gainloss_freq_all$expand == "Non Expanded","cellsize_pval"]=NA
gainloss_freq_all$cellsize_padj[1:88] = p.adjust(gainloss_freq_all$cellsize_pval[1:88],method="fdr") 
gainloss_freq_all$forcolor = factor(paste0(gainloss_freq_all$gain_loss,"-",gainloss_freq_all$expand))
exp_color= c("red","blue","grey","grey");
names(exp_color) = c("Gain-Expanded","Loss-Expanded","Gain-Non Expanded","Loss-Non Expanded")

gainloss_freq_all$label = -log10(gainloss_freq_all$binary_padj)
gainloss_freq_all$label = unlist(lapply(gainloss_freq$label,function(x){
  if(x >=3){ # p < 1e-3
    "***"
  }else if(x >=2){ # p < 1e-2
    "**"
  }else if(x >=-log10(0.05)){# p < 0.05
    "*"
  }else{
    NA
  }
}))
gainloss_freq$label_position = gainloss_freq$ratio + 0.08
pdf("cnv_clone_level/expanded_vs_nonexpanded_cnv_barplot.pdf",width=5,height=5)
bar <- ggplot(gainloss_freq_all,aes(y=chrarm,x=ratio,fill=forcolor)) + 
geom_bar(stat="identity",position="dodge") + facet_wrap(~gain_loss) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     text=element_text(size=10))+ 
  xlab("Ratio") + ylab("") + xlim(0, max(gainloss_freq_all$ratio)+0.15)+
  scale_fill_manual(values=exp_color) + 
  labs(fill="Gain/Loss") 
bar
dev.off()

saveRDS(gainloss_freq_all)

