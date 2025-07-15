#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=rem_regions
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

while read sample
do
	sample_path=$(echo $sample | sed "s/.*\/${TCGA_project}\(.*\)\/variants.vcf.gz/\1/g")
	sample_path="/g/strcombio/fsupek_cancer3/NMD/TCGA_germline_variants/"${TCGA_project}$sample_path
	if [ -f "${sample_path}/${TCGA_project}_regions_filt1.vcf.gz" ]; then
		echo "skipping -->" ${sample_path}
	else
		sbatch remove_regions_sample_vcf.sh $TCGA_project $sample
		#sleep 30
	fi
done < ~/data/TCGA_germline_variants/${TCGA_project}_VCFs_samples_list.txt

