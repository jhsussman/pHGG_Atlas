# Clonal Analysis Description

This folder contains the code that generates Fig.
6a and Extended Fig.9.

Overall, Clonalscope was used to estimate CNVs and detect subclones in the scRNA-Seq datasets.
For each patient, an early time point sample is chosen as the "seed", and Clonalscope utilizes the seed to lineage trace the subclones in all the rest samples belonging to the same patient.

**Summary of each script:**

1.  **smooth_snRNA_seg_inferCNV_normal.R**: script for running Clonalscope for each patient's each sample. InferCNV and manually checked normal cells from the integrative atlas is used as the reference/control cells for Clonalscope.

2.  **run_all_sample_smooth.sh**: bash script for submitting (1) on a server.

3.  **snRNA_lineage_trace_inferN_smooth_seg.R**: script for lineage tracing given (1)'s earliest timepoint result as seed/prior.

4.  **run_lineage_trace_inferN_smoothSeg.sh**: bash script for submitting (3) on a server.

5.  **plot_functions.R**: containing helper functions for converting CNVkit segmentation and creating bulk CNV level diagnostic plots.

6.  **full_match.txt**: matching between each sample and their corresponding WGS file.

7.  **lineage_match_progressiveSeed_14patient.txt:** metadata including patient id, selected earliest timepoint sample of each patient as the prior for lineage tracing, and all the rest sample name of each corresponding patient.

8.  **plot_help_functions.R**: containing code for calculating tumor percentages and creating fishplots for lineage tracing results.

9.  **Clonal_fishplot_SharedScripts.R *(for Fig. E9a)*:** create fish plots of tumor clonal percentages at different timepoints for each patient. Define expanded clones.

10. **Clonal_DEGs_SharedScripts.R:** Perform DEG on each tumor clone versus rest tumor cells in one patient; Summarize & plot top DEGs among all patients. Also calculates gene set/pathway enrichment analysis of expanded clones based on the clonal DEG.

11. **plot_scRNA_CNV_summary.R *(for Fig.6a, Fig. E9b,c)*:** Summarizes Clonalscope CNV profile & lineage tracing results at chromosomal arm level for subclones across all patients. Calculated subclone CNV correlation and also the CNV ratio bar plots.
