#!/usr/bin/bash

DIR="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"
VCF_type=$1

while read cancer
do
echo $cancer;

commands_file="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_${VCF_type}_PTC_NMD_rules_samples_commands.sh"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

if [ $array_jobs -gt 1000 ]; then
	array_jobs_half=$((array_jobs/2 +1))
	split -l $array_jobs_half ${commands_file} /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_${VCF_type}_PTC_NMD_rules_samples_commands_
	sbatch --array=0-${array_jobs_half}:5%40 ${DIR}/NMD_rules_and_efficiency_array.sh /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_${VCF_type}_PTC_NMD_rules_samples_commands_aa
	sbatch --array=0-${array_jobs_half}:5%40 ${DIR}/NMD_rules_and_efficiency_array.sh /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_${VCF_type}_PTC_NMD_rules_samples_commands_ab 
else 
	sbatch --array=0-${array_jobs}:5%40 ${DIR}/NMD_rules_and_efficiency_array.sh ${commands_file}
fi
done < /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt
