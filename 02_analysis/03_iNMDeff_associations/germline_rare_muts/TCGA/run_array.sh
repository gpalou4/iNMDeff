#!/usr/bin/bash

commands_file="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/TCGA/sbatch_commands.txt"
commands_file_tmp="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/rare_germline_variants_associations/TCGA/sbatch_commands_"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

# Split file in subfiles of 1000 lines
if  [ $array_jobs -gt 999 ]; then
	split -l 999 ${commands_file} ${commands_file_tmp}
	for file in $(ls $commands_file_tmp*); 
	do 
        echo $file
		array_jobs_split=$(wc -l $file | awk '{print $1}') 
		sbatch --array=0-${array_jobs_split}%60 slurm_array.sh $file
	done
else 
	sbatch --array=0-${array_jobs}%60 slurm_array.sh ${commands_file}
fi
