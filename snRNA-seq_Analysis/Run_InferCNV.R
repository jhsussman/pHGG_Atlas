#!/usr/bin/env Rscript

options(error = function() traceback(2))

library("infercnv")

patientID <- c(
"C15498",
"C1027419",
"C1037505",
"C1060383",
"C1061121",
"C107625",
"C176874",
"C2399853",
"C2542041",
"C2542533",
"C2751264",
"C34809",
"C547104",
"C70848",
"C714384",
"C799746")

for(i in c(6,15)) {
# create the infercnv object
  sample <- patientID[i]
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=sprintf("./matrix_count_%s.txt.gz", sample),
                                    annotations_file=sprintf("./cell_table_%s.txt", sample),
                                    delim="\t",
                                    gene_order_file="./gencode_downsampled.v3.txt",
                                    ref_group_names=c("White Blood Cells","Vascular Cells")
                                  )

  out_dir=sprintf("./output_dir_%s", sample)
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=FALSE, 
                             plot_steps=FALSE,
                             analysis_mode="subclusters",
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads=60
                             )

}
