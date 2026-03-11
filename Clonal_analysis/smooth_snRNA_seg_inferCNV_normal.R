library(Clonalscope)
library(EnsDb.Hsapiens.v86)
library(Seurat)
library(ggplot2)
library(stringr)

setwd("~/nzhanglab/data/CPTCA_Derek_pHGG/Derek_snRNA_20221219/") # set path to the clonalscope folder
source("~/nzhanglab/project/jrong/CPTCA_pHGG/analysis/plot_functions.R")

analysis_folder <-"/home/stat/jrong/nzhanglab/project/jrong/CPTCA_pHGG/"
#args=c("HTAN_pHGG_3403_Reg1_snRNA","HTAN_pHGG_3403_Reg2")
#args=c("HGG_7316_2295_R1","7316_2295_reg1") # sample 2
args=commandArgs(trailingOnly = TRUE)
start.time <- Sys.time()
sample = args[2] #"HTAN_pHGG_371_Reg2" #"7316_1057_reg2"
snRNA_name = args[1]#"HTAN_pHGG_371_Reg2_snRNA" #"HGG_7316_1057_R2"
dir_path = paste0(analysis_folder,"inferN_smooth_results/smooth_seg_tumorOnly2/",snRNA_name)  # output or results folder
if((is.na(sample)) | (nchar(sample) <1)){
  paired_bulk=FALSE
}else{
  paired_bulk=TRUE
}
filter_WGS_cor = FALSE
if(filter_WGS_cor){
  dir_path=paste0(analysis_folder,"smooth_seg_results/",snRNA_name,"/filter_WGS_cor/")
}
if (!file.exists(dir_path)){dir.create(file.path(dir_path),recursive=T)}
seg_prop=0.05;ngene_filter=100;#ngene_filter=150


# load the integrated snRNA-Seq data
snrna_obj <- readRDS("snRNA_merged_SeuratObj_QCFiltered_DoubletFiltered_downsampled_rPCA_integrated_clustered_withInferCNVTumorNormal_withNeftelTypes.rds")
# Size of each chromosome (hg19 and GRCh38 are provided.)
size=read.table("/home/stat/jrong/nzhanglab/project/jrong/clonal_scope/Clonalscope/data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)
# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS("/home/stat/jrong/nzhanglab/project/jrong/clonal_scope/Clonalscope/data-raw/cyclegenes.rds")
# bed file indicating gene positions (hg19 and GRCh38 are provided.)
bed=read.table("/home/stat/jrong/nzhanglab/project/jrong/clonal_scope/Clonalscope/data-raw/hg38.genes.bed", sep='\t', header = T)
gtf=readRDS("/home/stat/nzh/nzhanglab/project/cywu/ref/refdata-gex-GRCh38-2020-A/genes/genes_filtered.rds")
gene_id=sapply(strsplit(sapply(strsplit(sapply(strsplit(gtf[,9],";"),'[',1),' '),'[',2),'[.]'),'[',1)
gtf$gene_id=gene_id
gene_name=sapply(strsplit(sapply(strsplit(gtf[,9],"gene_name "),'[',2),';'),'[',1)
bed=gtf[,c(1,4,5,10)]

# create CNV segments from WGS
if(paired_bulk){
  cnvkit_cns = read.table(paste0("CNVkit/",sample,".call.cns"),header=T,sep="\t") 
  # iteratively smoothing the segments
  smoothed_cns = smooth_cns_segment(cnvkit_cns,seg_prop=seg_prop)
  #layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
  #layout(mat = layout.matrix,
  #       heights = c(5, 5), # Heights of the two rows
  #       widths = c(10, 10)) # Widths of the two columns
  #plotbulkCNV(cnvkit_cns,title="Original WGS Segments")
  #plotbulkCNV(smoothed_cns,title=paste0("Smoothed Segments with prop=",seg_prop))
  
  # Generate segmentation table for each chromosome arm (w/ paired WGS)
  seg_table_filtered=data.frame("chr"=str_remove(smoothed_cns$chromosome,"chr"),
                                'start'=smoothed_cns$start,
                                'end'=smoothed_cns$end,
                                'states'=smoothed_cns$cn/2, # copy number state
                                'length'=smoothed_cns$end-smoothed_cns$start,
                                'mean'=0, 'var'=0, 'Var1'=1:nrow(smoothed_cns),'Freq'=50000,stringsAsFactors = F)
  seg_table_filtered$chrr=paste0(seg_table_filtered[,1],":", seg_table_filtered[,2])
}else{ # no paired bulk, using chromosomal arms as a large segments
  chrarm=read.table("/home/stat/jrong/nzhanglab/project/jrong/clonal_scope/Clonalscope/data-raw/cytoarm_table_hg38.txt", stringsAsFactors = F, sep='\t', header=T)
  chrarm=chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
  bin_bed=chrarm[,-2]
  seg_table_filtered=data.frame("chr"=bin_bed[,1], 'start'=as.numeric(bin_bed[,2]),
                                'end'=as.numeric(bin_bed[,3]), 'states'=1, 'length'=as.numeric(bin_bed[,3])-as.numeric(bin_bed[,2]),
                                'mean'=0, 'var'=0, 'Var1'=1:nrow(bin_bed),'Freq'=50000,
                                'chrr'=paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)
}

# extract corresponding sample matrix info
mtx = snrna_obj@assays$RNA@counts[,snrna_obj$orig.ident == snRNA_name] # previously using RNA
#mtx = snrna_obj@assays$integrated@counts[,snrna_obj$orig.ident == snRNA_name] 
barcodes = as.data.frame(colnames(mtx))
features=data.frame(rownames(mtx),stringsAsFactors = F)
features=cbind(gene_id[match(features[,1], gene_name)], features[,1])

# remove cycle genes
Input_filtered <- FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)
# create cell barcode with matched (optional, set )
# celltype = data.frame(barcode = colnames(mtx),
#                      celltype=snrna_obj@active.ident[snrna_obj$orig.ident == snRNA_name])
# rownames(celltype)=celltype[,1] 
# Remove raw inputs
rm(mtx); rm(barcodes); rm(features);

# load final tumor annotation
metadata=snrna_obj@meta.data
rm(snrna_obj)
tumor_anno <- readRDS("../Wenbao_scATAC/seurat_rna_Neuroglia_final.rds")
# load inferCNV normal-tumor annotation for Clonalscope
celltype = data.frame(barcode = colnames(Input_filtered$mtx),
                      celltype=metadata[colnames(Input_filtered$mtx),"TumorNormal"])
# using final tumor annotation
rownames(celltype)=celltype[,1] 
celltype[celltype[,2] == "Normal",2] = "normal"
celltype[celltype[,2] == "Tumor",2] = "tumor"
celltype$TumorAnno = metadata[colnames(Input_filtered$mtx),"TumorNormal"]
tumor_subclone_bc=intersect(rownames(celltype),rownames(tumor_anno@meta.data[tumor_anno$orig.ident %in% snRNA_name,]))
celltype[tumor_subclone_bc,"TumorAnno"] =  tumor_anno$cell_state[tumor_subclone_bc]

rm(tumor_anno)
include='tumor'; # want tumor only
mtx=Input_filtered$mtx;barcodes=Input_filtered$barcodes;features=Input_filtered$features;
celltype0=celltype; var_pt=0.99; var_pt_ctrl=0.99; 
alpha_source='all'; ctrl_region=NULL; 
seg_table_filtered=seg_table_filtered;size=size;
dir_path=dir_path; breaks=50; prep_mode = 'intersect'; est_cap = 2;
alpha=2; beta=2; niter=200; sigmas0=NULL; U0=NULL; Z0=NULL;
clust_mode='all';
clustering_barcodes=NULL; seed=200; clustering0=NULL; result0=NULL; verbose=FALSE;
burnin=NULL; thinning=1; mincell = NULL; cutoff = 0.2;
threshold_2nd=-0.2; re_est=NULL; Est_read1=FALSE; Est_read2=FALSE; Clust_read1=FALSE; Clust_read2=FALSE;

plot_path=paste0(dir_path,"/cov_hist.pdf")
plot_path2=paste0(dir_path,"/cov_hist_updated.pdf")
result_all=list()
# parameter check and reset
if(is.null(re_est)){
  if(is.null(celltype0)){
    re_est=TRUE
  }else if(!is.null(celltype0)){
    re_est=FALSE
  }
}
if(is.null(clustering_barcodes)){
  clustering_barcodes=barcodes[,1]
}
if(is.null(celltype0)){
  celltype0=cbind(barcodes[,1], rep("normal", nrow(barcodes)))
}
# estimate deltas (change of read coverage) in each segment
if( Est_read1==TRUE){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
  deltas_all=readRDS(paste0(dir_path, "/deltas_allseg.rds"))
}else{
  deltas_all=EstRegionCov(mtx=mtx, barcodes=barcodes, features=features, bed=bed, 
                          celltype0=celltype0, var_pt=var_pt, var_pt_ctrl=var_pt_ctrl, 
                          ngene_filter = ngene_filter,include=include,
                          alpha_source=alpha_source, ctrl_region=ctrl_region, 
                          seg_table_filtered=seg_table_filtered, size=size,
                          plot_path=plot_path, breaks=breaks)
  
  saveRDS(deltas_all, paste0(dir_path, "/deltas_allseg.rds"))
}

# Filter segments by the direction of WGS CNV state and gene expression
if(filter_WGS_cor){
  #dir_path=paste0("debug_results/",snRNA_name,"/filter_WGS_cor/")
  cell_seg_cov_change = deltas_all$deltas_all
  a = log(unlist(lapply(cell_seg_cov_change,mean)))
  keep_seg_idx = ((a * (deltas_all$seg_table_filtered$states - 1)) >= 0)
  kept_cns = seg_table_filtered[keep_seg_idx,]
  colnames(kept_cns)[4] = "cn"
  
  deltas_all_filtered = list()
  deltas_all_filtered$deltas_all = deltas_all$deltas_all[keep_seg_idx]
  deltas_all_filtered$ngenes = deltas_all$ngenes[keep_seg_idx]
  deltas_all_filtered$seg_table_filtered = deltas_all$seg_table_filtered[keep_seg_idx,]
  deltas_all_filtered$alpha = deltas_all$alpha
  deltas_all_filtered$barcodes_normal = deltas_all$barcodes_normal
  df_obj=PrepCovMatrix(deltas_all=deltas_all_filtered, ngene_filter=50, prep_mode=prep_mode)
  #plotbulkCNV(smoothed_cns,title="WGS Filtering by Corr_Direction",filtered_seg=kept_cns,filter_color="#0000ff50")
}else{
  df_obj=PrepCovMatrix(deltas_all=deltas_all, ngene_filter=ngene_filter, prep_mode=prep_mode)
  # plot and save the segments being selected
  if(paired_bulk){
    pdf(paste0(dir_path,"/bulk_WGS.pdf"))
    layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
    layout(mat = layout.matrix,
         heights = c(5, 5), # Heights of the two rows
         widths = c(10, 10)) # Widths of the two columns
    plotbulkCNV(cnvkit_cns,title="Original WGS Segments")
    plotbulkCNV(smoothed_cns,filtered_seg=df_obj$seg_table_filtered,title="Smoothed and Filtered Segments")
    dev.off()
  }
}
df=df_obj$df
seg_table_filtered=df_obj$seg_table_filtered
print(paste0("Dimension after filtering:",dim(df)))

cna_states_WGS=seg_table_filtered$states
cna_states_WGS[which(cna_states_WGS>1)]= 1.5
cna_states_WGS[which(cna_states_WGS<1)]= 0.5
cna_states_WGS=as.numeric(cna_states_WGS)
df2=apply(df, c(2), function(x) pmin(x, est_cap))

# With prior CNV segment results, estimate cluster
if(!is.null(clustering0)){
  if(is.null(result0)){
    clustering02=MCMCtrim(clustering0)
    result0=AssignCluster(clustering02)
  }else{
    if(ncol(clustering0$data)==ncol(result0$Uest)){
      result0=result0
    }else{
      stop("clustering0 and result0 do not match!")
    }
  }
  ## subset the regions
  df2=df2[,which(colnames(df2) %in% colnames(clustering0$data)), drop=F]
  U0=result0$Uest
  U0=U0[,which(colnames(clustering0$data) %in% colnames(df2)), drop=F]
  cna_states_WGS=clustering0$priors$cna_states_WGS[which(colnames(clustering0$data) %in% colnames(df2))]
}

# Non-parametric clustering step
if( Clust_read1==TRUE){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
  clustering=readRDS(paste0(dir_path, "/nonpara_clustering.rds"))
}else{
  if(clust_mode=='all'){
    selr=1:length(cna_states_WGS)
    # clustering=BayesNonparCluster(Xir=df2[which(rownames(df2) %in% clustering_barcodes), , drop=F], cna_states_WGS =cna_states_WGS , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0[which(rownames(df2) %in% clustering_barcodes)] , seed = seed)
  }else if(clust_mode=='cna_only'){
    selr=which(as.character(cna_states_WGS)!='1')
  }else{
    stop("Please specify a valid clust_mode.")
  }
  clustering=BayesNonparCluster(Xir=df2[which(rownames(df2) %in% clustering_barcodes),selr], cna_states_WGS =cna_states_WGS[selr] , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0[which(rownames(df2) %in% clustering_barcodes)] , seed = seed, verbose=verbose)
  
  saveRDS(clustering,paste0(dir_path, "/nonpara_clustering.rds"))
}

# cut clusters
clustering2=MCMCtrim(clustering, burnin = burnin, thinning = thinning)

## sel mincell
if(is.null(mincell)){
  mincell=round(length(clustering_barcodes)*0.01)
}
# merge clusters
result=AssignCluster(clustering2, mincell = mincell, cutoff = cutoff)
Zest=result$Zest
print(paste0("Zest: ",table(Zest)))

df_obj$df=df_obj$df[which(rownames(df_obj$df) %in% clustering_barcodes),, drop=F]

result1=list( deltas_all= deltas_all, cna_states_WGS=cna_states_WGS, df_obj=df_obj, clustering=clustering, clustering2=clustering2, result=result)
result_all[['result1']]=result1

# update coverage estimation if needed
result2=NULL
tmp=1
while((any(result$corrs< threshold_2nd) & tmp<10) & re_est==TRUE){
  tmp=tmp+1
  message("Estimated clusters with negative correlation to the WGS data.")
  message(paste0("Start estimation ",tmp))
  
  new_normal=rownames(clustering2$data)[which(result$corrs %in% unique(result$corrs)[which(unique(result$corrs)< threshold_2nd)])]
  celltype0[which(celltype0[,2]=='normal'),2]='normal_pre'
  celltype0[which(celltype0[,1] %in% new_normal),2]='normal'
  
  if( Est_read2==TRUE & paste0("deltas_allseg_updated",tmp,".rds") %in% list.files(dir_path)){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
    deltas_all=readRDS(paste0(dir_path, "/deltas_allseg_updated",tmp,".rds"))
  }else{
    deltas_all=EstRegionCov(mtx=mtx, barcodes=barcodes, features=features, bed=bed, celltype0=celltype0, var_pt=var_pt, var_pt_ctrl=var_pt_ctrl, include=include,
                            alpha_source=alpha_source, ctrl_region=ctrl_region, seg_table_filtered=seg_table_filtered, size=size,
                            plot_path=plot_path2, breaks=breaks)
    
    saveRDS(deltas_all, paste0(dir_path, "/deltas_allseg_updated",tmp,".rds"))
  }
  df_obj=PrepCovMatrix(deltas_all=deltas_all, ngene_filter=ngene_filter)
  df=df_obj$df
  print(paste0("Dimension after filtering:",dim(df)))
  seg_table_filtered=df_obj$seg_table_filtered
  
  cna_states_WGS=seg_table_filtered$states
  cna_states_WGS[which(cna_states_WGS>1)]= 1.5
  cna_states_WGS[which(cna_states_WGS<1)]= 0.5
  cna_states_WGS=as.numeric(cna_states_WGS)
  
  df2=apply(df, c(2), function(x) pmin(x, est_cap))
  
  if(!is.null(clustering0)){
    # clustering02=MCMCtrim(clustering0)
    # result0=AssignCluster(clustering02)
    
    ## subset the regions
    df2=df2[,which(colnames(df2) %in% colnames(clustering0$data)), drop=F]
    U0=result0$Uest
    U0=U0[,which(colnames(clustering0$data) %in% colnames(df2)), drop=F]
    cna_states_WGS=clustering0$priors$cna_states_WGS[which(colnames(clustering0$data) %in% colnames(df2))]
  }
  
  if( Clust_read2==TRUE & paste0("nonpara_clustering_updated",tmp,".rds") %in% list.files(dir_path)){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
    clustering=readRDS(paste0(dir_path, "/nonpara_clustering_updated",tmp,".rds"))
    
  }else{
    if(clust_mode=='all'){
      selr=1:length(cna_states_WGS)
      #clustering=BayesNonparCluster(Xir=df2, cna_states_WGS =cna_states_WGS , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0 , seed = seed)
    }else if(clust_mode=='cna_only'){
      selr=which(as.character(cna_states_WGS)!='1')
    }else{
      stop("Please specify a valid clust_mode.")
    }
    clustering=BayesNonparCluster(Xir=df2[which(rownames(df2) %in% clustering_barcodes),selr], cna_states_WGS =cna_states_WGS[selr] , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0[which(rownames(df2) %in% clustering_barcodes)] , seed = seed, verbose=verbose)
    
    saveRDS(clustering,paste0(dir_path, "/nonpara_clustering_updated",tmp,".rds"))
  }
  clustering2=MCMCtrim(clustering, burnin = burnin, thinning = thinning)
  
  result=AssignCluster(clustering2, mincell = mincell, cutoff = cutoff)
  table(result$Zest)
  Zest=result$Zest
  
  print(paste0("Zest: ",table(Zest)))
  
  df_obj$df=df_obj$df[which(rownames(df_obj$df) %in% clustering_barcodes),, drop=F]
  
  result2=list( deltas_all= deltas_all, cna_states_WGS=cna_states_WGS, df_obj=df_obj, clustering=clustering, clustering2=clustering2, result=result, new_normal=new_normal)
  result_all[[paste0('result',tmp)]]=result2
  
  if(tmp==10){
    warning("After 10 rounds of iterations, there is still a cluster with corrs<-0.3.")
  }
  #Est_read2=FALSE
  #Clust_read2=FALSE
}

# save results
result_all[['celltype0']]=celltype0
ll=sort(as.numeric(gsub("result","",names(result_all)[grepl("result", names(result_all))])))[length(names(result_all)[grepl("result", names(result_all))])]
result_all[['result_final']]=result_all[[paste0('result',ll )]]

Cov_obj = result_all

clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
table(result$Zest)
saveRDS(Cov_obj,paste0(dir_path,paste0("/","Cov_obj",".rds")))

# save the clustering result related plots
pdf(paste0(dir_path,paste0("/Clustering_Results",".pdf")))
if(paired_bulk){
  layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
  layout(mat = layout.matrix,
       heights = c(5, 5,5), # Heights of the two rows
       widths = c(10, 10,10)) # Widths of the two columns
  plotbulkCNV(cnvkit_cns,title="Original CNV Segments")
  plotbulkCNV(smoothed_cns,title="Smoothed CNV Segments")
 #plotbulkCNV(smoothed_cns,title="WGS Correlation Filtering",filtered_seg=kept_cns,filter_color="#0000ff50")
  plotbulkCNV(smoothed_cns,filtered_seg=df_obj$seg_table_filtered,title="Smoothed Segments with Selected (Red)")
}
g <- PlotClusters(df = Cov_obj$result_final$df_obj$df, celltype = celltype, #celltype0, 
                  Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode='annot')
print(g)
# Plot the UMAP and save the coordinates
emb=PlotUMAP(df = Cov_obj$result_final$df_obj$df, celltype = celltype, Assign_obj =result, mode = "Zest")
dev.off()


