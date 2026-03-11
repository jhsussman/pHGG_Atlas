#!/bin/bash
#$ -N smooth_seg
#$ -l m_mem_free=70G
#$ -o /home/stat/jrong/nzhanglab/project/jrong/CPTCA_pHGG/logs/individual_sample/smooth_seg/inferCNV_normal_tumorOnly 
#$ -j y

data_file=/home/stat/jrong/nzhanglab/data/CPTCA_Derek_pHGG/Derek_snRNA_20221219/full_match.txt

file_lists=$(cat $data_file)
if [[ -n $SGE_TASK_ID ]]; then
  F=$(cat $data_file | sed -n ${SGE_TASK_ID}p)
fi

eval "array=($F)"
arg1="${array[0]}" # snRNA name
arg2="${array[1]}" # bulk sample name

echo $arg1 
echo $arg2

Rscript ../smooth_snRNA_seg_inferCNV_normal.R $arg1 $arg2
