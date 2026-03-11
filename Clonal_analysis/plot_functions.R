## visualization

# create x axis index of CNV segments
makefull=function(x,nr,name){
  xnew= sum(nr)
  from=1
  for(i in 1:length(x)){
    to = from+nr[i]-1
    # assign CNV state
    xnew[from:to] = x[i]
    # assign chr segment length
    names(xnew)[from:to] = name[i]
    from=to+1
  }
  xnew
}

# adapted from Clonalscope https://github.com/seasoncloud/Clonalscope/blob/main/R/Segmentation_bulk.R
####
# @ cnvkit_cns: original CNVkit output
# @ title: plot title
# @ filtered_seg: segments filtered by Clonalscope's standard
####
plotbulkCNV <-function(cnvkit_cns,title=NULL,filtered_seg=NULL,filter_color="#FF000050"){
  #dev.new(width=width, height=height,unit=unit) 
  # code adpated from Nancy's explore_cancer_discovery_niches.R
  start=as.numeric(cnvkit_cns$start)
  end=as.numeric(cnvkit_cns$end)
  chrname=cnvkit_cns$chromosome
  chrname=strsplit(chrname, "chr")
  chrname=unlist(lapply(chrname, function(x){x[2]}))
  len=end-start
  nreps=ceiling(len/1e6)

  # plot segments based on segment length
  clustfull=makefull(cnvkit_cns$cn, nreps,paste0(chrname,":",cnvkit_cns$start))
  plot(clustfull, col="black",xaxt="n",xaxs="i", xlab="", ylab="",pch=20, cex=0.3,main=title)
  title(xlab="Chromosome", line=2, cex.lab=1.2)
  title(ylab="Copy Number States", line=2, cex.lab=1.2)
  abline(h=2, lwd=2, col="blue"); grid()
  xval0=0
  for(i in 1:length(nreps)){
    xval=sum(nreps[1:i])
    if(i==length(nreps) ||chrname[i+1] != chrname[i]) {
      segments(xval,0, xval, max(cnvkit_cns$cn), col="blue", lty=2)
      mtext(chrname[i], side=1, at=(xval0+xval)/2)
      xval0=xval
    }
  }
  
  # draw the filtered segments on the same plot
  if(!is.null(filtered_seg)){
    # find original index of the segments
    for (chrr in filtered_seg$chrr){
      idx = which(names(clustfull) == chrr)
      xleft=min(idx); xright=max(idx);
      rect(xleft,0,xright,max(cnvkit_cns$cn),col=filter_color,border = NA)
    }
  }
}


# This functions smooths the bulk CNV segments from CNVkit
smooth_cns_segment <- function(cnvkit_cns,seg_prop=NULL){
  cnvkit_cns$cn[cnvkit_cns$cn > 3] = 3
  cnvkit_cns$cn[cnvkit_cns$cn < 2] = 1
  smoothed_cns = c()
  start_idx = 1
  end_idx=1
  chr = "chr1"
  # merge segments with same copy states
  for(i in 1:dim(cnvkit_cns)[1]){
    #for(i in 1:550){
    #print(c(chr,start_idx,end_idx))
    next_chr = cnvkit_cns[i+1,1] # chromosome
    if(!is.na(cnvkit_cns[i+1,"cn"]) & (next_chr==chr) & (cnvkit_cns[i,"cn"]==cnvkit_cns[i+1,"cn"])){
      end_idx = end_idx + 1
    }else{ # switching chrom or ending of file
      smoothed_cns = rbind(smoothed_cns,c(chr,
                                          cnvkit_cns[start_idx,"start"],cnvkit_cns[end_idx,"end"],
                                          cnvkit_cns[end_idx,"cn"]))
      # update start and end idx
      start_idx = i+1
      end_idx = i+1
      chr = next_chr
    }
  }
  smoothed_cns=as.data.frame(smoothed_cns)
  colnames(smoothed_cns) = c("chromosome","start","end","cn")
  smoothed_cns[, colnames(smoothed_cns)[2:4]] <- lapply(colnames(smoothed_cns)[2:4], function(x) as.numeric(smoothed_cns[[x]]))
  
  # further combining segments based on neighbors and size
  if(!is.null(seg_prop)){
    merged_cns = smoothed_cns
    #merged_cns$seg_len = merged_cns$end - merged_cns$start
    temp_cns = c()
    # Iteratively merge & smooth every time
    while(TRUE){ 
      #print("iter")
      change = 0
      for(i in 2:(dim(merged_cns)[1]-1)){
        #if(i %in% skip_list){next}
        #for(i in 2:dim(merged_cns)[1]){
        chr = merged_cns[i,1] # current chromosome
        pre_chr = merged_cns[i-1,1] # previous chromosome
        next_chr = merged_cns[i+1,1] # next chromosome
        # current segment proportion in itself + nearby segments
        cur_prop = (merged_cns[i,3]-merged_cns[i,2])/(merged_cns[i,3]-merged_cns[i,2]+merged_cns[i-1,3]-merged_cns[i-1,2]+merged_cns[i+1,3]-merged_cns[i+1,2])
        # (1) merging segments if its neighbors have same cn state 
        # and if it is small compared to all 3 segments
        if(!is.na(merged_cns[i+1,"cn"]) & 
           (pre_chr==chr) & (next_chr==chr) & (merged_cns[i-1,"cn"]==merged_cns[i+1,"cn"]) &
           (cur_prop <= seg_prop)){
          
          temp_cns = rbind(merged_cns[0:(i-2),],
                           c(chr,merged_cns[i-1,2],merged_cns[i+1,3],merged_cns[i-1,"cn"]))
          colnames(temp_cns)= c("chromosome","start","end","cn")
          temp_cns = rbind(temp_cns,merged_cns[(i+1):dim(merged_cns)[1],])
          temp_cns=as.data.frame(temp_cns)
          colnames(temp_cns) = c("chromosome","start","end","cn")
          temp_cns[, colnames(temp_cns)[2:4]] <- lapply(colnames(temp_cns)[2:4], function(x) as.numeric(temp_cns[[x]]))
          #print(i)
          #print("Yes")
          change =change+1
          break
        }
      }
      # then smooth temporary data matrix
      merged_cns = c()
      start_idx = 1
      end_idx=1
      chr = "chr1"
      # merge segments with same copy states
      for(i in 1:dim(temp_cns)[1]){
        #for(i in 1:550){
        #print(c(chr,start_idx,end_idx))
        next_chr = temp_cns[i+1,1] # chromosome
        if(!is.na(temp_cns[i+1,"cn"]) & (next_chr==chr) & (temp_cns[i,"cn"]==temp_cns[i+1,"cn"])){
          end_idx = end_idx + 1
        }else{ # switching chrom or ending of file
          merged_cns = rbind(merged_cns,c(chr,temp_cns[start_idx,"start"],temp_cns[end_idx,"end"],
                                          temp_cns[end_idx,"cn"]))
          # update start and end idx
          start_idx = i+1
          end_idx = i+1
          chr = next_chr
        }
      }
      merged_cns=as.data.frame(merged_cns)
      colnames(merged_cns) = c("chromosome","start","end","cn")
      merged_cns[, colnames(merged_cns)[2:4]] <- lapply(colnames(merged_cns)[2:4], function(x) as.numeric(merged_cns[[x]]))

      # if no new segment change, stop merging
      if(change ==0){
        break
      }      
    }
  }
  # returning results
  if(!is.null(seg_prop)){
    return(merged_cns)
  }else{
    return(smoothed_cns)
  }
}

# This function divides the table based on say ~ every 5Mb segments
# and the identity is determined by majority vote + smoothing
weighted_vote_cns_segment <- function(cnvkit_cns,seg_prop=NULL,seg_size=5000000){
  cnvkit_cns$cn[cnvkit_cns$cn > 3] = 3
  cnvkit_cns$cn[cnvkit_cns$cn < 2] = 1
  # smoothing firsting
  smoothed_cns = smooth_cns_segment(cnvkit_cns)
  temp_cns = c()
  # for each chromosome
  for(chr in paste0("chr",c(1:22,"X","Y"))){
    # divide by segment length
    chr_idx = which(smoothed_cns$chromosome == chr)
    seg_len = (smoothed_cns$end - smoothed_cns$start)[chr_idx]
    #cum_len = cumsum(seg_len)
    cum_seg_len = 0 # cumulative segment length
    seg_list = c()
    for(idx in chr_idx){
      cum_seg_len = cum_seg_len + (smoothed_cns$end[idx] - smoothed_cns$start[idx])
      seg_list = append(seg_list,idx)
      # add in each segment that adds up >= segment size (e.g. 1Mb)
      if(cum_seg_len >= seg_size){
        # calculate weighted copy number state for each segment size region
        weighted_cn =round(sum(smoothed_cns$cn[seg_list]*(smoothed_cns$end[seg_list] - smoothed_cns$start[seg_list])/sum(smoothed_cns$end[seg_list] - smoothed_cns$start[seg_list])))
        new_seg = c(chr,smoothed_cns[seg_list[1],"start"],
                    smoothed_cns[tail(seg_list,1),"end"], weighted_cn)
        temp_cns = rbind(temp_cns,new_seg)
        cum_seg_len = 0
        seg_list = c()
      }else if(idx == tail(chr_idx,1)){
        weighted_cn =round(sum(smoothed_cns$cn[seg_list]*(smoothed_cns$end[seg_list] - smoothed_cns$start[seg_list])/sum(smoothed_cns$end[seg_list] - smoothed_cns$start[seg_list])))
        new_seg = c(chr,smoothed_cns[seg_list[1],"start"],
                    smoothed_cns[tail(seg_list,1),"end"], weighted_cn)
        temp_cns = rbind(temp_cns,new_seg)
      }
    }
  }
  temp_cns=as.data.frame(temp_cns)
  colnames(temp_cns) = c("chromosome","start","end","cn")
  temp_cns[, colnames(temp_cns)[2:4]] <- lapply(colnames(temp_cns)[2:4], function(x) as.numeric(temp_cns[[x]]))
  # smooth again after vote by segments
  smoothed_cns = smooth_cns_segment(temp_cns)
  return(smoothed_cns)
}

### This function plots out correlation plot of clonalscope and inferCNV results
plot_cor <- function(clonalscope_vec,inferCNV_vec){
  # calculate correlation
  # convert to one-hot encoding
  cor_df = data.frame(clonalscope=clonalscope_vec,inferCNV=inferCNV_vec)
  dmy <- dummyVars(" ~ .", data = cor_df)
  trsf <- data.frame(predict(dmy, newdata = cor_df))
  # create correlation heatmap
  mydata.rcorr = rcorr(as.matrix(trsf))
  mydata.cor = cor(trsf, method = c("pearson"))
  corrplot(mydata.cor,method = 'number', col=colorRampPalette(c("blue","white","red"))(200))
}


### This function calculates correlation of with CNV
calc_gene_cnv_corr <- function(){
  ## compute new means
  Usub=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
  for(kk in ind_clusters){
    meanX=colMeans(Xir[which(Zest==kk),, drop=F],  na.rm = T)
    meanX[is.na(meanX)]=0
    Usub[kk,]=meanX
  }
  # compute correlation
  corrs=apply(Usub2[,corr_region_ind,drop=F], 1, function(x){
  tmp=sort(((as.numeric(x)-1)*(as.numeric(priors[corr_region_ind])-1))/(sqrt(sum((as.numeric(x)-1)^2))*sqrt(sum((as.numeric(priors[corr_region_ind])-1)^2))))
  })
}

### This function calculates Tumor vs Normal counts given a segment
# code adapted from Clonalscope EstCovRegion.R
# return index of seg table that shall be kept
calc_seg_wilcoxon <- function(mtx,features,Cov_obj,ngene_filter=50){
  rownames(mtx)  = features[,1]
  seg_table_filtered = Cov_obj$result_final$df_obj$seg_table_filtered
  # annotate cells based on final annotation
  celltype0=data.frame(barcode=colnames(mtx),celltype="Filtered")
  rownames(celltype0) = colnames(mtx)
  celltype0[names(Cov_obj$result_final$result$Zest)[Cov_obj$result_final$result$annot == "T"],2] ="tumor"
  celltype0[names(Cov_obj$result_final$result$Zest)[Cov_obj$result_final$result$annot == "N"],2] ="normal"
    
  rna_var=apply(mtx,1, var)
  ngenes=Matrix::rowSums(mtx)
  rna=mtx[which(rna_var<quantile(rna_var,var_pt) & (rna_var!=0) & ngenes>ngene_filter), ]
  rna_control=mtx[which(rna_var<quantile(rna_var,var_pt_ctrl) & (rna_var!=0) & ngenes>ngene_filter), ]
  
  bed_sub=bed[match(rownames(rna), bed[,4]),, drop=F]
  bed_sub[is.na(bed_sub)]=0
  
  sel_ind=c()
  i = 1
  # calculate gene count for each segment
  for(chrr in paste0(seg_table_filtered$chrr)){
    query=GRanges(seqnames = paste0('chr',seg_table_filtered$chr[which(seg_table_filtered$chrr==chrr)]),
                  ranges = IRanges(as.numeric(seg_table_filtered$start[which(seg_table_filtered$chrr==chrr)]),
                                   as.numeric(seg_table_filtered$end[which(seg_table_filtered$chrr==chrr)])))
    subject=GRanges(seqnames = bed_sub[,1], ranges=IRanges(start=as.numeric(bed_sub[,2]), end=as.numeric(bed_sub[,3])))
    ov=as.matrix(findOverlaps(query=query, subject = subject))
    gene_ind=ov[,2]
    rna_sub=rna[gene_ind,, drop=F]
    sel_cell=which(Matrix::colSums(rna_sub)>0)
    
    if(length(sel_cell)==0){
      next
    }
    #sel_ind=c(sel_ind, chrr)
    rna_sub=rna_sub[,sel_cell, drop=F]
    
    celltypes=celltype0
    celltypes=celltypes[match(colnames(rna), celltypes[,1]),]
    celltypes=celltypes[sel_cell,, drop=F]
    # normal cell counts
    barcodes_normal=colnames(rna_sub)[which(celltypes[,2]=='normal')]
    controlCounts=as.matrix(rna_sub[,which(celltypes[,2]=='normal'), drop=F])
    # tumor cell counts
    #testCounts=as.matrix(rna_sub[,which(celltypes[,2]=='tumor'), drop=F]) 
    testCounts=as.matrix(rna_sub)
    
    # compute cell size from control region
    alphas=Matrix::colSums(rna_control[,sel_cell, drop=F])
    alpha_controls=alphas[which(celltypes[,2]=='normal')]/median(alphas)  #mean
    alpha_tests=alphas/median(alphas)
    
    testCounts=testCounts/matrix(rep(alpha_tests, nrow(testCounts)), nrow=nrow(testCounts), byrow = T)
    controlCounts=controlCounts/matrix(rep(alpha_controls, nrow(controlCounts)), nrow=nrow(controlCounts), byrow = T)
    deltas=Matrix::colSums(testCounts)/median(Matrix::colSums(controlCounts))
    
    normal_deltas = deltas[which(celltypes[,2]=='normal')]
    tumor_deltas = deltas[which(celltypes[,2]=='tumor')]
    # wilcoxon rank sum test based on WGS CNV state
    cnv_state = seg_table_filtered$states[i]
    if(cnv_state > 1){
      pval = (wilcox.test(tumor_deltas, normal_deltas, alternative = "greater"))$p.value
      if(pval <= 0.025){sel_ind = append(sel_ind,i)}
    }else if(cnv_state < 1){
      pval = (wilcox.test(tumor_deltas, normal_deltas, alternative = "less"))$p.value
      if(pval < 0.025){sel_ind = append(sel_ind,i)}
    }else{ # copy neutral
      pval = (wilcox.test(tumor_deltas, normal_deltas, alternative = "two.sided"))$p.value
      if(pval < 0.05){sel_ind = append(sel_ind,i)}
    }
    i=i+1
  }
  return(sel_ind)
}


