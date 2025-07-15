#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=TCGA_array_commands_SBATCH
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

while read cancer
do
echo $cancer;

commands_file="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_endogenous_sbatch_commands_samples.txt"
commands_file_tmp="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_endogenous_sbatch_commands_samples_final_"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

# Split file in subfiles of 1000 lines
if  [ $array_jobs -gt 999 ]; then
	split -l 990 ${commands_file} ${commands_file_tmp}
	for file in $(ls $commands_file_tmp*); 
	do 
		array_jobs_split=$(wc -l $file | awk '{print $1}') 
		sbatch --array=0-${array_jobs_split}:20%10 endogenous_NMD_efficiency_regression_array.sh $file
	done

else 
	sbatch --array=0-${array_jobs}:20%10 endogenous_NMD_efficiency_regression_array.sh ${commands_file}
fi
done < /g/strcombio/fsupek_home/gpalou/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt
#done < TCGA_cancer.txt
