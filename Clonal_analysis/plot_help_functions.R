# This function plots out fishtail & clevRis plots 
plot_clonal_expansion <- function(i,lineage_match,metadata,snrna_obj,tumor_anno,percent=99.5,plot_path=NULL,plot=T,return=F,vlabSize=4){
  if(is.null(plot_path)){
    plot_path="plots/lineages/interesting_patients/Rest_patients/"
  }
  patient_id = lineage_match$Patient.ID[i]
  print(paste0(i,", patientID: ",patient_id ))
  seed_snRNA = lineage_match$snRNA_seed[i]
  snRNA_lineages = lineage_match$snRNA_fileAccessions[i]
  snRNAs = c(seed_snRNA,unlist(str_split(snRNA_lineages,",")))
  timepoints_all = metadata[snRNAs,c("timepoint")]
  plot_df_celltype= as.data.frame(cbind(barcodes = names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs]),
                                      snRNA_accession = snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs],
                                      celltype=snrna_obj$cellType[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]))
  plot_df_celltype$Clonalscope_cluster="Filtered/Unknown"
  plot_df_celltype$TumorAnno = snrna_obj$TumorNormal[names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs])]
  tumor_subclone_bc = intersect(rownames(tumor_anno),
                              names(snrna_obj$old.ident[snrna_obj$old.ident %in% snRNAs]))# then tumor subclones
  plot_df_celltype[tumor_subclone_bc,"TumorAnno"] = tumor_anno[tumor_subclone_bc,"cell_state"]

  Cov_obj_list <- list()
  # load Clonalscope result
  for (j in 1:length(snRNAs)){
  # load Clonalscope clustering result
    f = snRNAs[j]
    if(j==1){
    # initial seed sample
      Cov_obj = readRDS(paste0("inferN_smooth_results/smooth_seg_tumorOnly2/",f,"/Cov_obj.rds"))
      Cov0_clusters = levels(as.factor(as.numeric(Cov_obj$result_final$result$Zest)))
      plot_df_celltype[names(Cov_obj$result_final$result$Zest),"Clonalscope_cluster"] = Cov_obj$result_final$result$Zest
    }else{
      Cov_obj = readRDS(paste0("./lineage_results_inferN_smooth_TumorOnly_allpatients//",patient_id,"/",f,"/Cov_obj.rds"))
      plot_df_celltype[names(Cov_obj$result_final$result$Zest),"Clonalscope_cluster"] = Cov_obj$result_final$result$Zest
    #clust_anno = Cov_obj$result_final$result$Zest
    #for(clu in unique(clust_anno)){
    #  if(!(clu %in% Cov0_clusters)){
    #    clust_anno[clust_anno == clu]=paste0(clu,".",j)
    #  }
    #}
      plot_df_celltype[names(Cov_obj$result_final$result$Zest),"Clonalscope_cluster"] = Cov_obj$result_final$result$Zest#Cov_obj$result_final$result$Zest
      Cov_obj_list[[f]] <- Cov_obj
    }
  }
  plot_df_celltype = plot_df_celltype[tumor_subclone_bc,]
  #plot_df_celltype$Clonalscope_cluster[is.na(plot_df_celltype$Clonalscope_cluster)] = "Filtered"

  # intial and progressive, two time points for C1027419
  timepoints <- c(1:length(unique(timepoints_all))) 
  # generate fraction table
  clones = sort(as.numeric(setdiff(levels(as.factor(plot_df_celltype$Clonalscope_cluster)),c("Filtered/Unknown"))))
  fracTable = matrix(rep(0,length(timepoints)*length(clones)),ncol=length(timepoints))
  for(i in 1:length(unique(timepoints_all))){
    tp = unique(timepoints_all)[i]
    tp_scRNA_name = snRNAs[timepoints_all==tp]
    for(j in 1:length(clones)){
      clone=clones[j]
      fracTable[j,i] = sum((plot_df_celltype$snRNA_accession %in% tp_scRNA_name) &
                           (plot_df_celltype$Clonalscope_cluster == clone))
    }
  }
  # normalize each column (timepoint) to sum up to 100
  fracTable = t(t(fracTable)/colSums(fracTable)) * percent # not * 100 to avoid bug in "createSeaObject"
  rownames(fracTable) <- clones
  colnames(fracTable) <- unique(timepoints_all)
  parents=rep(0,length(clones))

  ##seaObject with enabled time point estimation
  seaObject_tp <- createSeaObject(fracTable, parents, timepoints,
                                timepointInterpolation = T,
                                cloneLabels=as.character(clones),
                                #col = rainbow(7)
  )
  # Dolphin plot
  if(plot==T){
    dir.create(paste0(plot_path,patient_id))
    pdf(paste0(plot_path,patient_id,"/Clone_expansion_dolphinplot_clevRvis.pdf"),
      height=5,width=10)
    par(mar=c(10,10,10,10))
    g<- dolphinPlot(seaObject_tp, showLegend = TRUE, main = patient_id, 
            vlines = timepoints, 
            vlab = unique(timepoints_all), #c("Initial CNS","Progressive(Non-Autopsy)","Autopsy"), 
            vlabSize = 4, 
            ylab = 'Cancer cell fraction', #annotations = annotsTable, 
            #pos = 'bottom', 
            separateIndependentClones = TRUE)
    print(g)
    dev.off()
  
    # Fishplot 
    pdf(paste0(plot_path,patient_id,"/Clone_expansion_fishplot.pdf"),
      height=7,width=14)
    #par(mar=c(2,2,2,2))
    #create a fish object
    fish = createFishObject(fracTable,parents,timepoints=timepoints,
                        clone.labels=as.character(clones),clone.annots=as.character(clones))
    # set up plotting colors for additional clones
    default_fishplot_cols = col = rev(c("#00008F", "#0000FF", "#0070FF", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#EF0000", "#888888"))
    if(length(clones) > length(default_fishplot_cols)){
      cols = colorRampPalette(brewer.pal(9, "Set1"))(length(clones))
      #setCol(fish,cols)
      fish = createFishObject(fracTable,parents,timepoints=timepoints,
                            clone.labels=as.character(clones),clone.annots=as.character(clones),col=cols)
    }
    #calculate the layout of the drawing
    fish = layoutClones(fish,separate.independent.clones=T)
    #draw the plot, using the splining method (recommended)
    #and providing both timepoints to label and a plot title
    fishPlot(fish,shape="spline",#title.btm=patient_id,
         cex.title=3, cex.vlab=1,
         vlines=timepoints,vlab=unique(timepoints_all),#c("Initial CNS","Progressive(Non-Autopsy)","Autopsy"),
         title=patient_id)

    dev.off()
  }else{
    par(mar=c(10,10,10,10))
    g<- dolphinPlot(seaObject_tp, showLegend = TRUE, main = patient_id, 
                    vlines = timepoints, 
                    vlab = unique(timepoints_all), #c("Initial CNS","Progressive(Non-Autopsy)","Autopsy"), 
                    vlabSize = vlabSize, 
                    ylab = 'Cancer cell fraction', #annotations = annotsTable, 
                    #pos = 'bottom', 
                    separateIndependentClones = TRUE)
  }
  if(return==T){
    return(g)
  }
}