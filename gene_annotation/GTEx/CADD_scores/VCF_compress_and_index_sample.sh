#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_annovar
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=2GB
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

GTEx_sample=$1
GTEx_samples_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs"

# Decompress
gunzip ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf.gz
# Compress with bgzip
conda activate samtools
bgzip ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf
# Index with tabix
tabix ${GTEx_samples_path}/${GTEx_sample}_liftOver_annotated.HG_multianno_rare_variants.vcf.gz

