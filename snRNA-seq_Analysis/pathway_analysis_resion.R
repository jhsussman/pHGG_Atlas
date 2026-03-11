library(data.table)
library(Seurat)
library(fgsea)
library(enrichR)
library(ggplot2)
library(viridis)

#degs = readRDS('data/intermediate/rna/degs_seurat_rna_Neuroglia_final_noIntegration_regrMITO.rds')
degs = readRDS('data/intermediate/rna/wilcox_degs_seurat_rna_Neuroglia_revision.rds')
degs = degs[avg_log2FC > 0.5 & pct.1 > pct.2]

## use all degs ####
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
degs$score = degs$avg_log2FC
degs = degs[order(-score)]
sele.dbs <- c("GO_Biological_Process_2023", 'KEGG_2021_Human',
              "Reactome_2022", 'MSigDB_Hallmark_2020')

cell_state = 'HSP+ Cells'

degs0 = degs[cluster == cell_state]

genes0 = unique(degs0$gene)

message(paste0(cell_state, ': ', length(genes0), ' genes'))
enriched_res = enrichr(genes0, sele.dbs)

if(is.null(enriched_res)){
  message('No enriched terms found!')
}else{
  for(db0 in names(enriched_res)){
    res.state0 = data.table(enriched_res[[db0]])
    res.state0[, 'Count' := as.numeric(unlist(strsplit(Overlap, '/', fixed = T))[1]), 
               by = Overlap]
    res.state0 = res.state0[order(P.value)]
    
    res.state0[, 'score' := -log10(P.value)]
    tmp = res.state0[Count >= 3 & P.value < 0.05]
    if(nrow(tmp) == 0) {
      stop(paste0(cell_state, ': no sig terms!')) 
      next
    }
    tmp = tmp[order(P.value)]
    nmax = min(15, nrow(tmp))
    pic.height = 8
    if(nmax < 5) pic.height = 2
    if(nmax >= 5 & nmax < 10) pic.height = 5
    
    
    tmp = tmp[1:nmax, ]
    tmp$Term = factor(tmp$Term, levels = rev(tmp$Term))
    p0 <- ggplot(tmp, aes(Term, score, fill = score)) +
      geom_bar(stat = 'identity') +
      coord_flip() + scale_fill_viridis(option = 'D') +
      labs(x="", title = paste0(db0, ': ', cell_state),
           y = '-log10(pvalue)') + guides(fill=guide_legend(title="-log10(pvalue)")) +
      theme_classic() + theme(axis.text.y = element_text(size = 12))
    #p0 <- plotEnrich(enriched_res, showTerms = 20, numChar = 40, y = "Count", 
    #                 orderBy = "P.value")
    cell_state1 = gsub(' ', '', cell_state)
    cell_state1 = gsub( '/', '-', cell_state1, fixed = T)
    cell_state1 = gsub( '+', '_', cell_state1, fixed = T)
    dir.create(paste0('Figures/RNA/enrichR/AllDEGs/', db0),
               showWarnings = F, recursive = T)
    filekey = paste0('Figures/RNA/enrichR/AllDEGs/', db0, '/', 
                     db0, '_', cell_state1)
    ggsave(p0, filename = paste0(filekey, '.pdf'),
           device = 'pdf', width = 12, height = pic.height)
    fwrite(res.state0[Count >= 3 & P.value < 0.05], sep = '\t',
           file = paste0(filekey, '.tsv'))
  }
  
}


