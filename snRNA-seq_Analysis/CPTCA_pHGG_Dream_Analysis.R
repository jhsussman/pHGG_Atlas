#### FINAL DREAM ANALYSIS
library('variancePartition')
library('edgeR')
library('BiocParallel')

# Read in starting Seurat object; either 1) "Tumor" for neoplastic cells only or 2) "Myeloid" for myeloid cells only
analysis_mode <- "Tumor"
if(analysis_mode == "Tumor") {
  so <- readRDS("/mnt/isilon/tan_lab/yuw1/R_work_dir/pHGG/SeuratObjects/seurat_rna_Neuroglia_final.rds")
  seq_depth_threshold <- 1e6
  cpm_threshold <- 10
  out_string <- ""
}
if(analysis_mode == "Myeloid") {
  so <- readRDS("/mnt/isilon/tan_lab/sussmanj/pHGG/snRNA-seq/Macrophages/Reanalysis_scRNA_macrophage_filtered2.rds")
  seq_depth_threshold <- 2e5
  cpm_threshold <- 50
  out_string <- "_myeloid"
}

# Correct autopsy annotation for sample 7316-6910
so$timepoint[so$orig.ident %in% c("HTAN_pHGG_6910_R1_scRNA", "HTAN_pHGG_6910_R2_scRNA")] <- "Progressive (Autopsy)"

# Clean timepoint coding
so$timepoint_n <- 2
so$timepoint_n[so$timepoint=="Initial CNS Tumor"] <- 1
so$timepoint_n[so$sampleID=="6477"] <- 3
so$timepoint_n[so$sampleID=="6151"] <- 1
so$timepoint_n[so$sampleID=="6513"] <- 2
so$timepoint_n <- factor(so$timepoint_n, levels=1:3)
table(so$patient_id, so$timepoint_n)

# Compile pseudobulk matrices by each unique sample
uniqSample <- sort(unique(so$orig.ident))
lPseudo <- list()
for(i in 1:length(uniqSample)) {
  print(i)
  idx <- so$orig.ident == uniqSample[i]
  lPseudo[[i]] <- apply(so@assays$RNA@counts[,idx,drop=F], 1, sum)
}
counts <- data.frame(do.call(cbind, lPseudo))
colnames(counts) <- uniqSample

# Compile sample annotation data frame
dm <- data.frame(
  patient   = tapply(so$patient_id, so$orig.ident, unique),
  timepoint = paste("T", tapply(so$timepoint_n, so$orig.ident, unique), sep=""),
  stringsAsFactors = F
)
dm$PatientTimepoint <- paste(dm$patient, dm$timepoint, sep=":")
mapSampleToTimepoint <- tapply(as.character(so$timepoint), so$orig.ident, unique)
dm$autopsy <- mapSampleToTimepoint[rownames(dm)] == "Progressive (Autopsy)"
rownames(dm) <- colnames(counts)

## Run Dream Analysis, ref: https://mirrors.dotsrc.org/bioconductor-releases/3.8/bioc/vignettes/variancePartition/inst/doc/dream.html

# Sample filtering; patient C15498 removed because of prior history of HGG before our "timepoint1"; patients C70848 and C547104 removed because missing paired data after QC filtering
colToKeep <- apply(counts, 2, sum) >= seq_depth_threshold 
colToKeep <- colToKeep & dm$timepoint %in% c("T1", "T2") 
colToKeep <- colToKeep & dm$patient != "C70848" & dm$patient != "C15498" & dm$patient != "C547104"
countsFiltered <- counts[,colToKeep]

# Filter lowly expressed genes using a CPM threshold of corresponding to 10 reads at minimum read count in at least 10 samples; ref: https://combine-australia.github.io/RNAseq-R/slides/RNASeq_filtering_qc.pdf
isexpr = rowSums(cpm(countsFiltered) > cpm_threshold) >= 10

geneExpr = DGEList( countsFiltered[isexpr,] )
geneExpr = calcNormFactors( geneExpr ) 

set.seed(8174651)

dmFiltered <- dm[colToKeep,]
design = model.matrix( ~ timepoint, dmFiltered)
vobj_tmp = voom( geneExpr, design, plot=TRUE)

param = SnowParam(30, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ timepoint + (1|patient) 

# Estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, dmFiltered, BPPARAM=param)

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, dmFiltered )
fitmm = eBayes(fitmm)
ttdream <- topTable( fitmm, coef='timepointT2', number=Inf)
rnk <- data.frame(gene=rownames(ttdream), t = ttdream$t, stringsAsFactors = F)
write.table(rnk,   sprintf("./CPTCA_pHGG_dream_mixed_effect_model_univariate_timepoint2%s_2024-09-10.rnk", out_string), quote=F, sep="\t", row.names=F, col.names = F)
write.csv(ttdream, sprintf("./CPTCA_pHGG_dream_mixed_effect_model_univariate_timepoint2%s_2024-09-10.csv", out_string), quote=F, row.names=T)

# Rerun the above Limma-Voom + Dream analysis, but adjust for autopsy effect
form <- ~ timepoint + autopsy + (1|patient) 

# Estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights( geneExpr, form, dmFiltered, BPPARAM=param )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, dmFiltered )
fitmm = eBayes(fitmm)

# Save autopsy-adjusted output
ttdream <- topTable( fitmm, coef='timepointT2', number=Inf)
rnk <- data.frame(gene=rownames(ttdream), t = ttdream$t, stringsAsFactors = F)
write.table(rnk,   sprintf("./CPTCA_pHGG_dream_mixed_effect_model_autopsyAdjusted_timepoint2%s_2024-09-10.rnk", out_string), quote=F, sep="\t", row.names=F, col.names = F)
write.csv(ttdream, sprintf("./CPTCA_pHGG_dream_mixed_effect_model_autopsyAdjusted_timepoint2%s_2024-09-10.csv", out_string), quote=F, row.names=T)

ttdream <- topTable( fitmm, coef='autopsyTRUE', number=Inf)
rnk <- data.frame(gene=rownames(ttdream), t = ttdream$t, stringsAsFactors = F)
write.table(rnk,   sprintf("./CPTCA_pHGG_dream_mixed_effect_model_autopsyAdjusted_autopsyTRUE%s_2024-09-10.rnk", out_string), quote=F, sep="\t", row.names=F, col.names = F)
write.csv(ttdream, sprintf("./CPTCA_pHGG_dream_mixed_effect_model_autopsyAdjusted_autopsyTRUE%s_2024-09-10.csv", out_string), quote=F, row.names=T)
