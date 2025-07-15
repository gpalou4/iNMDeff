#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=rem_regions_VCF
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=12GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

# Activate base environment
conda activate pandas_bcftools

TCGA_project=$1
sample=$2

# Create sample directories if do not exist
sample_path=$(echo $sample | sed "s/.*\/${TCGA_project}\(.*\)\/variants.vcf.gz/\1/g")
sample_path="/g/strcombio/fsupek_cancer3/NMD/TCGA_germline_variants/"${TCGA_project}$sample_path

if [ ! -f "${sample_path}/${TCGA_project}_regions_filt1.vcf.gz" ]; then
	mkdir -p $sample_path
	# Remove regions
	vcftools --gzvcf $sample --exclude-bed FASTA_BED.ALL_GRCh38.novel_CUPs.bed --recode --recode-INFO-all --stdout | gzip -c > $sample_path/${TCGA_project}_regions_filt1.vcf.gz
fi
