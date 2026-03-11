#!/bin/bash
#$ -N tumor_lineage
#$ -l m_mem_free=120G
#$ -o /home/stat/jrong/nzhanglab/project/jrong/CPTCA_pHGG/logs/lineage/interesting
#$ -j y

data_file=/home/stat/jrong/nzhanglab/data/CPTCA_Derek_pHGG/Derek_snRNA_20221219/lineage_match_progressiveSeed_allFinal.txt

file_lists=$(cat $data_file)
if [[ -n $SGE_TASK_ID ]]; then
  F=$(cat $data_file | sed -n ${SGE_TASK_ID}p)
fi

eval "array=($F)"
arg1="${array[0]}" # snRNA name
arg2="${array[1]}" # bulk sample name
arg3="${array[2]}"
arg4="${array[3]}"

echo $arg1 
echo $arg2
echo $arg3
echo $arg4

# Rscript snRNA_lineage_trace.R $arg1 $arg2 $arg3 $arg4
Rscript ../snRNA_lineage_trace_inferN_smooth_seg.R $arg1 $arg2 $arg3 $arg4
