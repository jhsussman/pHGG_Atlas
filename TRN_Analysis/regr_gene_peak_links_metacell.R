
source('scDataAnalysis_Utilities_simp.R')
`%notin%` = Negate(`%in%`)
library(parallel)
library(Signac)

## get a binary matrix indicates the gene-peak affinity
## gene.list are data.table including column name gene_name 
## gene_ann should include gene_name,chr,start,end
get_gene2peak_map <- function(gene.list, peak_names,
                              gene_ann, distal_dist = 2e05){
  
  peaks <- tidyr::separate(data.table(peak_name = peak_names), 
                           col = peak_name, remove = F,
                           into = c('chr', 'start', 'end'))
  class(peaks$start) = class(peaks$end) = 'integer'
  setkey(gene_ann, gene_name)
  setkey(peaks, peak_name)
  gene.list = gene.list[gene_name %in% gene_ann$gene_name]
  gene.list$chr = gene_ann[gene.list$gene_name]$chr
  gene.list$start = gene_ann[gene.list$gene_name]$start
  gene.list$end = gene_ann[gene.list$gene_name]$end
  
  ## for each gene, get the corresponding peaks
  gene2peaks = lapply(gene.list$gene_name, function(x) {
    
    chr0 = gene_ann[x]$chr
    start0 = gene_ann[x]$start
    end0 = gene_ann[x]$end
    
    peaks0 = peaks[chr == chr0]
    peaks0 = peaks0[abs(start/2 + end/2 - start0/2 - end0/2) <= distal_dist]
    return(peaks0$peak_name)
  } )
  
  ## pool all peaks relate to one gene
  gene2peaks.u <- lapply(sort(unique(gene.list$gene_name)), function(x){
    id = which(gene.list$gene_name == x)
    tmp_peak <- do.call('c', lapply(id, function(x) gene2peaks[[x]]))
    return(tmp_peak)
  })
  names(gene2peaks.u) <- sort(unique(gene.list$gene_name))
  lens = sapply(gene2peaks.u, length)
  
  genes.f <- names(which(lens > 0))
  lens = lens[lens > 0]
  ## construct overlap matrix
  gene2peaks.dt <- data.table('gene' = rep(genes.f, lens),
                              'peak' = do.call('c', lapply(genes.f, 
                                                           function(x) gene2peaks.u[[x]])))
  upeaks = sort(unique(gene2peaks.dt$peak))
  gene2peaks.dt[, 'id1' := which(genes.f == gene), by = gene]
  gene2peaks.dt[, 'id2' := which(upeaks == peak), by = peak]
  gene2peak.map <- sparseMatrix(i = gene2peaks.dt$id1,
                                j = gene2peaks.dt$id2,
                                dimnames = list(genes.f, upeaks))
  gene2peak.map = gene2peak.map * 1
  
  return(gene2peak.map)
}

## annotate peaks with gene +/- 5kb of its TSS
# input peak_coords with chr-start-end, format
annPeak2Gene <- function(peak_coords, gene_ann, proxim_dist = 5e+03){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))
  peaks = data.table(peaks)
  
  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    
    peaks0$gene_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + proxim_dist) & tss >= (peaks0$start[i] - proxim_dist)]
      if(nrow(tss0) > 0 ) {
        peaks0$gene_name[i] = paste(unique(tss0$gene_name), collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  setkey(peaks_ann, peak_name)
  
  return(peaks_ann)
  
}


## map gene to overlapping atac peak
## gene_list with genename, chr, start, end
geneOverlapPeak <- function(gene_list, peak_names, mid_dist = 1000){
  # should include tss information in gene_list
  peaks = tidyr::separate(data = data.table('peak_name' = peak_names),
                          col = peak_name, into = c('chr', 'start', 'end'),
                          remove = F)
  peaks = data.table(peaks)
  class(peaks$chr) = 'character'
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    gene0$peak_name = 'Not_Found'
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id0' := any( abs(peaks0$midP -start) < mid_dist | abs(peaks0$midP - end) < mid_dist), 
          by = gene_name]
    gene1 = gene0[peak_id0 == FALSE]
    gene2 = gene0[peak_id0 == TRUE]
    gene2[, 'peak_id' :=  which.min(abs(peaks0$midP - start - 1000)), by = gene_name]
    
    gene2[, 'peak_name' :=  peaks0[peak_id]$peak_name, by = gene_name]
    gene2$peak_id = NULL
    gene_new = rbind(gene_new, gene1, gene2)
  }
  gene_new[, c('peak_id0') := NULL]
  return(gene_new)
}


## preparing inputs just need to do one time
filterMetacell_Done = FALSE # if true, skip combining and filtering step
if(!filterMetacell_Done){
  ## combine all coembed results from a single patient ####
  seuratAtacPath = 'SeuratObjects/seurat_atac_signacIntegrated_bySample_withReference_4036_4037.rds'
  seurat.atac = readRDS(seuratAtacPath)
  sampleNames = unique(seurat.atac$sampleID)
  doCombine = T # just run this once
  metacell.rna.mtx = metacell.atac.mtx <- list()
  if(doCombine){
    for(sampleName in sampleNames){
      tmp1 <- readRDS(paste0('Coembed_Results/metacells/rna_metacell_mtx_',
                             sampleName, '.rds'))
      tmp2 <- readRDS(paste0('Coembed_Results/metacells/atac_metacell_mtx_',
                             sampleName, '.rds'))
      metacell.rna.mtx[[sampleName]] = tmp1
      metacell.atac.mtx[[sampleName]] = tmp2
    }
  }
  
  metacell.rna.mtx = do.call('cbind', metacell.rna.mtx)
  metacell.atac.mtx = do.call('cbind', metacell.atac.mtx)
  dim(metacell.rna.mtx)
  dim(metacell.atac.mtx)
  
  ## filter the features of metacells ####
  ## filter peaks that accessible in less than 5% of all cell type
  peaks.mean.ctype <- sapply(unique(seurat.atac$ATAC_snn_res.0.4), function(x){
    rowMeans(seurat.atac@assays$ATAC@counts[, seurat.atac$ATAC_snn_res.0.4 == x] > 0)
  })
  
  rmax = apply(peaks.mean.ctype, 1,  max)
  summary(rmax)
  filtered.peaks = names(which(rmax > 0.02))
  
  ## further filtered peaks by dap or fc peaks
  dap_bulk = readRDS(file = paste0('data/intermediate/all_atac_pseudo_bulk_limma_voom.rds'))
  fc_bulk = readRDS('data/intermediate/atac_pseudo_bulk_byTumorState_topFC.rds')
  open.frac.pks = readRDS('MetaData/peak_openFrac_byState.rds')
  
  filtered.peaks = intersect(filtered.peaks, c(dap_bulk$peak, fc_bulk$peak,
                                                   rownames(open.frac.pks)))
  cand.peaks = rownames(metacell.atac.mtx)
  names(cand.peaks) = NULL
  
  all(filtered.peaks %in% cand.peaks)
  metacell.atac.mtx = metacell.atac.mtx[cand.peaks %in% filtered.peaks, ]
  
  ## focus on selected genes (or sele.genes)
  sele.genes = intersect(rownames(metacell.rna.mtx), rownames(seurat.atac[['ACTIVITY']]))
  # filter genes
  seurat.rna <- readRDS('SeuratObjects/seurat_rna_chop.rds')
  mtx = seurat.rna@assays$RNA@counts
  mtx.ctype = sapply(unique(seurat.rna$cellType), function(x) rowMeans(mtx[, seurat.rna$cellType == x] >0) )
  mfreq_gene = rowMax(mtx.ctype)
  sele.genes = intersect(sele.genes, rownames(mtx)[mfreq_gene > 0.05])
  
  ## gene-peak-affinity map
  ## construct gene peak affinity binary matrix
  gene_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/gene_ann_hg38.txt')
  gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]
  
  gene2peak.map <- get_gene2peak_map(gene.list = data.table('gene_name' = sele.genes),
                                     peak_names = filtered.peaks, 
                                     gene_ann = gene_ann,
                                     distal_dist = 5e5)
  
  peaks.used = colnames(gene2peak.map)
  sele.genes = sele.genes[sele.genes %in% rownames(gene2peak.map)]
  
  
  
  metacell.rna.mtx = metacell.rna.mtx[sele.genes, ]
  metacell.atac.mtx = metacell.atac.mtx[peaks.used, ]
  
  save(metacell.rna.mtx, metacell.atac.mtx, gene2peak.map,
       file = 'Coembed_Results/metacells/rna_atac_metacell_matrices.RData')
  
  
}

## do regression gene by gene ####
load('Coembed_Results/metacells/rna_atac_metacell_matrices.RData')
sele.genes = rownames(metacell.rna.mtx)
peaks.used = rownames(metacell.atac.mtx)

metacells.mtx = rbind(metacell.rna.mtx, metacell.atac.mtx)
metacells.df = data.frame(t(metacells.mtx))
pnames = names(metacells.df)[-c(1:length(sele.genes))]
pnames = gsub('.', '-', pnames, fixed = T)
names(metacells.df)[c(1:length(sele.genes))] = sele.genes
names(metacells.df)[-c(1:length(sele.genes))] = pnames

## do it in parallel
runParallel = FALSE 
stime = Sys.time()

if(runParallel){
  regress_function <- function(gene0) {
    peaks = names(which(gene2peak.map[gene0, ] == 1))
    data <- subset(metacells.df, select = c(gene0, peaks))
    names(data)[1] = 'y'
    model <- lm(y ~ -1 + ., data=data)
    return(coef(summary(model)))
  }
  # Use the mclapply function to perform the regressions in parallel
  ncores = detectCores() - 1
  regr_list <- mclapply(sele.genes,  FUN = regress_function, 
                        mc.cores = ncores)
  message(paste('#cores:', ncores))
  names(regr_list) = sele.genes
  
}else{
  curate_dataList = list()
  for(gene0 in sele.genes){
    peaks = names(which(gene2peak.map[gene0, ] == 1))
    data <- subset(metacells.df, select = c(gene0, peaks))
    names(data)[1] = 'y'
    curate_dataList[[gene0]] = data
  }
  
  regr_list = list()
  for(x in sort(sele.genes)){
    rdata <- curate_dataList[[x]]
    regr_list[[x]] <- coef(summary(lm(y ~ . -1, data = rdata)))
  }
}
etime = Sys.time()
etime-stime
saveRDS(regr_list, "EP_Prediction/regrRes4ep_prediction_metacell.rds")

## summarize/filter loops ####
regr_list = readRDS("EP_Prediction/regrRes4ep_prediction_metacell.rds")
regr.sum <- lapply(names(regr_list), function(t){
  x = regr_list[[t]]
  x = x[, c(1, 4), drop = F]
  x = data.frame(x)
  x = data.table(x, keep.rownames = T)
  x$gene_name = t
  return(x)
})

regr.sum = do.call('rbind', regr.sum)
names(regr.sum)[c(1, 3)] = c('peak_name', 'P_value')
regr.sum[, 'p_val_bonf' := pmin(1, P_value*nrow(regr.sum))]
regr.sum$fdr = p.adjust(regr.sum$P_value, method = 'fdr')

regr.sum$peak_name = sapply(regr.sum$peak_name, function(x) gsub('`', '', x, fixed = T)  )
predicted_links = regr.sum[fdr < 0.05 & abs(Estimate) > 0.1 & grepl(peak_name, pattern = '^chr')]

fwrite(predicted_links, file = 'EP_Prediction/regrRes4_gene_peak_links_metacell.tsv',
       sep = '\t')

useEnhancerPeakOnly = T
if(useEnhancerPeakOnly){
  ##older version
  if(F){
    gene_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/gene_ann_hg38.txt')
    gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]
    
    gene_ann_sele = gene_ann[gene_name %in% predicted_links$gene_name, ]
    
   
    ## avoid peaks to overlap with promoters -- this step is optional 
    #-- all transcripts
    tss_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/transcript_ann_hg38.txt')
    tss_ann = tss_ann[gene_biotype %in% c('protein_coding', 'lincRNA', 'miRNA')]
    tss_ann[, 'Tss' := ifelse(strand == '+', start, end)]
    
  }else{
    
    #--new version, use longest transcripts
    annotation = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/scGEL/data/annotation_hg38.rds')
    tss_ann <- Signac:::CollapseToLongestTranscript(annotation)
    tss_ann <- as.data.frame(resize(tss_ann, width = 1, fix = 'start'))
    names(tss_ann)[1] = 'chr'
    tss_ann$Tss = tss_ann$start
    gene_ann = tss_ann = data.table(tss_ann)
    
    gene_ann_sele = gene_ann[gene_name %in% predicted_links$gene_name, ]
    
  }
  
  setkey(gene_ann_sele, gene_name)
  predicted_links[, 'chr' := gene_ann_sele[J(predicted_links$gene_name)]$chr]
  predicted_links[, 'Tss' := gene_ann_sele[J(predicted_links$gene_name)]$Tss]
  predicted_links[, 'promoter_pos' := paste(chr, Tss - 1000, Tss + 1000, sep = '-')]
  
  # any peak within promoter region
  peak.ann = annPeak2Gene(unique(predicted_links$peak_name), tss_ann, 1000)
  setkey(peak.ann, peak_name)
  peaks.nprom = peak.ann[nchar(gene_name) == 0]$peak_name
  peaks.prom = peak.ann[nchar(gene_name) > 0]$peak_name
  
  predicted_links.ep = predicted_links[peak_name %in% peaks.nprom]
  
  ## assign nearest peak to promoter (within 1kb)
  gene_list <- subset(predicted_links.ep, select = c(gene_name, chr, Tss)) %>%
    .[!duplicated(.)]
  gene_list[, 'start' := Tss - 1000]
  gene_list[, 'end' := Tss + 1000]
  gene_list2peak = geneOverlapPeak(gene_list, peak_names = peaks.used,
                                   mid_dist = 1000)
  setkey(gene_list2peak, gene_name)
  predicted_links.ep[, 'promoter_peak' := gene_list2peak[J(predicted_links.ep$gene_name)]$peak_name]
  
  predicted_links.ep = predicted_links.ep[promoter_peak != 'Not_Found'] # filter genes whose promters are not open
  predicted_links.ep = subset(predicted_links.ep, select = c(gene_name, promoter_pos, 
                                                             promoter_peak, peak_name, Tss,
                                                             P_value, p_val_bonf, fdr, Estimate))
    
  predicted_links.ep[, 'start' := as.integer(unlist(strsplit(peak_name, '-'))[2]), by = peak_name]
  predicted_links.ep[, 'end' := as.integer(unlist(strsplit(peak_name, '-'))[3]), by = peak_name]
  
  predicted_links.ep[, 'ep_dist' := abs(Tss - start/2 - end/2)]
  predicted_links.ep = subset(predicted_links.ep, select = c(gene_name, promoter_pos, 
                                                             promoter_peak, peak_name, Tss, ep_dist, 
                                                             P_value, p_val_bonf, fdr, Estimate))
  #names(predicted_links.ep)[4] = 'enhancer_peak'
  fwrite(predicted_links.ep, file = 'EP_Prediction/regrRes4_gene_peak_ep_links_metacell.tsv',
         sep = '\t')
}
