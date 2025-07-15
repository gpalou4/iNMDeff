#!/usr/bin/bash

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx"
VCF_type=$1

rm /g/strcombio/fsupek_home/gpalou/out_err_files/NMDrules_array-*

while read GTEx_tissue
do
echo $GTEx_tissue;

commands_file="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/${GTEx_tissue}/${GTEx_tissue}_${VCF_type}_PTC_NMD_rules_samples_commands.sh"
commands_file_tmp="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/germline_PTCs/tissues/${GTEx_tissue}/${GTEx_tissue}_${VCF_type}_PTC_NMD_rules_samples_commands_"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

if [ $array_jobs -gt 1000 ]; then
	array_jobs_half=$((array_jobs/2 +1))
	split -l $array_jobs_half ${commands_file} ${commands_file_tmp}
	sbatch --array=0-${array_jobs_half}:10%20 ${DIR}/NMD_rules_and_efficiency_array.sh ${commands_file_tmp}_aa
	sbatch --array=0-${array_jobs_half}:10%20 ${DIR}/NMD_rules_and_efficiency_array.sh ${commands_file_tmp}_ab 
else 
	sbatch --array=0-${array_jobs}:10%20 ${DIR}/NMD_rules_and_efficiency_array.sh ${commands_file}
fi
done < /g/strcombio/fsupek_cancer1/gpalou/GTEx/RNAseq/GTEx_tissues_clean.txt
