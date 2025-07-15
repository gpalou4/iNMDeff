#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=rem_SNPs_%
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=6GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

# Activate base environment
conda activate pandas_bcftools

TCGA_project=$1
counter=0

sample_path=$(echo $sample | sed "s/.*\/${TCGA_project}\(.*\)\/variants.vcf.gz/\1/g")
sample_path="/g/strcombio/fsupek_cancer3/NMD/TCGA_germline_variants/"${TCGA_project}$sample_path

while read sample
do
	num_variants_after=$(zcat $sample_path/${TCGA_project}_regions_filt1.vcf.gz | grep -v "#" | wc -l)
	num_variants_before=$(zcat $sample | grep -v "#" | wc -l)
	echo | awk -v var1="$num_variants_before" -v var2="$num_variants_after" '{print (1-(var2/var1))*100}'

done < ~/data/TCGA_germline_variants/${TCGA_project}_VCFs_samples_list.txt
