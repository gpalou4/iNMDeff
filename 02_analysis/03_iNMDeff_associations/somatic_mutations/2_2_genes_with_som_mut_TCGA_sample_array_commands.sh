#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=array_TCGA_associations
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=4GB
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

TCGA_cancer=$1

OUT_DIR="/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/${TCGA_cancer:5}"
# Commands file
commands_file="${OUT_DIR}/${TCGA_cancer}_genes_som_mut_array_commands.txt"
commands_file_tmp="${OUT_DIR}/${TCGA_cancer}_genes_som_mut_array_commands_"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

# Split file in subfiles of 1000 lines
if  [ $array_jobs -gt 999 ]; then
	split -l 999 ${commands_file} ${commands_file_tmp}
	for file in $(ls $commands_file_tmp*); 
	do 
		array_jobs_split=$(wc -l $file | awk '{print $1}') 
		sbatch --array=0-${array_jobs_split}:15%30 2_genes_with_som_mut_TCGA_sample_array.sh $file
	done
else 
	sbatch --array=0-${array_jobs}:15%30 2_genes_with_som_mut_TCGA_sample_array.sh ${commands_file}
fi
