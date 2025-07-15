while read cancer
do
echo $cancer;

commands_file="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_ASE_samples_commands.sh"
array_jobs=$(wc -l $commands_file | awk '{print $1}')

if [ $array_jobs -gt 1000 ]; then
	array_jobs_half=$((array_jobs/2 +1))
	split -l $array_jobs_half ${commands_file} /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_ASE_samples_commands_
	sbatch --array=1-${array_jobs_half}%10 ASE_calling_sample_array.sh /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_ASE_samples_commands_aa
	sbatch --array=1-${array_jobs_half}%10 ASE_calling_sample_array.sh /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${cancer}/${cancer}_ASE_samples_commands_ab 
else 
	sbatch --array=1-${array_jobs}%10 ASE_calling_sample_array.sh ${commands_file}
fi
done < ~/data/TCGA_RNAseq_quantification/TCGA_projects_names.txt