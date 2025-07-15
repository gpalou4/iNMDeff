#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_liftover_sample
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=10GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment

GTEx_sample=$1

# Paths and files
GTEx_VCFs_path="/g/strcombio/fsupek_cancer1/gpalou/GTEx/WES_VCF/VCFs"
bad_regions="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/liftover/regions_to_remove/FASTA_BED.ALL_GRCh37.novel_CUPs.bed"
GTEx_sample_VCF=${GTEx_VCFs_path}/"${GTEx_sample}.vcf.gz"
GTEx_sample_VCF_corrected=${GTEx_VCFs_path}/"${GTEx_sample}_chr.vcf.gz"
GTEx_sample_VCF_corrected_final=${GTEx_VCFs_path}/"${GTEx_sample}_chr_badRegionsFilt.vcf.gz"
GTEx_sample_liftover_VCF=${GTEx_VCFs_path}/"${GTEx_sample}_liftOver.vcf.gz"
GTEx_sample_liftover_rejected_VCF=${GTEx_VCFs_path}/"${GTEx_sample}_liftOver_rejectedVariants.vcf.gz"

# 1) Add "chr" to the VCF
zcat ${GTEx_sample_VCF} | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | gzip -c > ${GTEx_sample_VCF_corrected}

# 2) Remove bad liftover regions
conda activate pandas_bcftools
vcftools --gzvcf ${GTEx_sample_VCF_corrected} --exclude-bed ${bad_regions} --recode --recode-INFO-all --stdout | gzip -c > ${GTEx_sample_VCF_corrected_final}

# 3) Perform liftover from hg19 to hg38
conda activate base

#if [ ! -f "${sample_path}/${TCGA_project}_regions_liftover_filt2.vcf.gz" ]; 
#then
java -Xmx10048m -jar /g/strcombio/fsupek_home/gpalou/software/picard/build/libs/picard.jar LiftoverVcf \
I=${GTEx_sample_VCF_corrected_final} \
O=${GTEx_sample_liftover_VCF} \
CHAIN=/g/strcombio/fsupek_home/gpalou/data/liftover/hg19ToHg38.over.chain \
REJECT=${GTEx_sample_liftover_rejected_VCF} \
R=/g/strcombio/fsupek_cancer1/gpalou/human_genome/GRCh38.p10.genome.fa \
WARN_ON_MISSING_CONTIG=true
#fi

# 4) Remove temporary files
rm ${GTEx_sample_VCF_corrected} ${GTEx_sample_VCF_corrected_final} 
