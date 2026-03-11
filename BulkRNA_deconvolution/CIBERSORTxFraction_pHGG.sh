#!/bin/bash
#SBATCH -c 4
#SBATCH --mem-per-cpu=80G
#SBATCH -t 48:00:00
#SBATCH --output=./batch_out/CIBERSORTxFraction_pHGG.out

## Note: fractons_latest.sif, which is the CIBERSORTx Fractions singularity image,
## must be downloaded to use with singularity command below

module load singularity

username=cibersort_username # must be registered to use CIBERSORTx
token=cibersort_token # must be registered to use CIBERSORTx
datadir=/mnt/isilon/tan_lab/tumultyj/MetaPathwayPaper/Brain/Glioma/pHGG/Deconvolution
ref=snRNA_cellState1CellTypeCondLabeled_counts.txt
mixture=bulkRNA_OpenPedCan_allPBTA_TPM.txt
savePre=OpenPedCan_allPBTA_cellState1CellTypeCondRef

singularity exec --bind /cm/shared --bind $datadir/data:/src/data --bind $datadir/outdir:/src/outdir/ --bind `pwd` fractions_latest.sif /src/CIBERSORTxFractions --username $username --token $token --single_cell TRUE --refsample $ref --mixture $mixture --rmbatchSmode TRUE --verbose TRUE --perm 500

cp $datadir/outdir/CIBERSORTx_Adjusted.txt $datadir/$savePre\_CIBERSORTx_Adjusted.txt
