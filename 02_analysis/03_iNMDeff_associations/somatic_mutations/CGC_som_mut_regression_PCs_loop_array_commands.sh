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
#SBATCH --time=15:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

NMD_method=$1 # Can be endogenous/PTC/ASE
type=$2 # "cancers" / "pancancer"
VAF=$3 # only for ASE
alpha=$4 # can be one of the alphas tested
robust_SPCA=$5 # yes/no

if [[ ${NMD_method} == "ASE" ]]; then
    NMD_method_char=${NMD_method}_${VAF}
else
    NMD_method_char=${NMD_method}
fi

commands_file="/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/array_files/${type}_${NMD_method_char}_CGC_som_mut_array_commands_${alpha}_robust_${robust_SPCA}.txt"
commands_file_tmp="/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/array_files/${type}_${NMD_method_char}_CGC_som_mut_array_commands_${alpha}_robust_${robust_SPCA}_"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

# Split file in subfiles of 1000 lines
if  [ $array_jobs -gt 999 ]; then
	split -l 999 ${commands_file} ${commands_file_tmp}
	for file in $(ls $commands_file_tmp*); 
	do 
		array_jobs_split=$(wc -l $file | awk '{print $1}') 
		sbatch --array=0-${array_jobs_split}:15%30 4_CGC_som_mut_regression_array.sh $file
	done
else 
	sbatch --array=0-${array_jobs}:15%30 4_CGC_som_mut_regression_array.sh ${commands_file}
fi
