
#!/usr/bin/bash

STAR_bam_sorted_file=$1
TCGA_BARCODE=$2
TCGA_PROJECT=$3
THREADS=$4
fastq_var=$5

#STAR_bam_sorted_file=$(readlink $STAR_bam_sorted_link)

# 1) Mapping quality: samtools -q
# ~10m
conda activate samtools

samtools view -h -b -q 255 $STAR_bam_sorted_file > ${TCGA_BARCODE}_q1.bam
samtools index ${TCGA_BARCODE}_q1.bam

# 2) Base quality:
# Strelka doesn't allow to use this filter, but internally it controls for this with their cutoff:
# In the germline SNV model, basecalls with mapping-adjusted qualities of 17 or less are filtered out. No such filtration is used by the somatic SNV model.

# 3) WASP method to correct for allele mapping bias
WASP_DIR="/g/strcombio/fsupek_home/gpalou/software/WASP/mapping"

# 3.1) Prepare input
# <2m

mkdir vcf

# VCF germline variants
sleep 30
VCF=$(grep ${TCGA_BARCODE} /local_scratch/gpalou/tmp/TCGA_germline_variants/${TCGA_PROJECT:5}_VCFs_samples_list.txt)
sleep 30
gunzip -d -c $VCF > vcf/germline_variants.vcf
bgzip -c vcf/germline_variants.vcf > vcf/germline_variants.vcf.gz
tabix -p vcf vcf/germline_variants.vcf.gz

# Split VCF per chromosome
conda activate bcftools

for i in {1..22}
do
        bcftools filter vcf/germline_variants.vcf.gz -r chr${i} -o vcf/chr${i}.vcf.gz -O z
done

for i in {X,Y,M}
do
        bcftools filter vcf/germline_variants.vcf.gz -r chr${i} -o vcf/chr${i}.vcf.gz -O z
done

# 3.2) get SNPs from VCF files:
# <1m

mkdir genotypes

${WASP_DIR}/extract_vcf_snps.sh vcf genotypes

# 3.3) Pull out reads that need to be remapped to check for allele mapping bias
# ~1h

conda activate /home/gpalou/anaconda3_envs/general
mkdir find_intersecting_snps

if [ ${fastq_var} == "single_end" ]; 
then
        python3 ${WASP_DIR}/find_intersecting_snps.py \
        --is_sorted \
        --output_dir find_intersecting_snps \
        --snp_dir genotypes \
        ${TCGA_BARCODE}_q1.bam
else
        python3 ${WASP_DIR}/find_intersecting_snps.py \
        --is_paired_end \
        --is_sorted \
        --output_dir find_intersecting_snps \
        --snp_dir genotypes \
        ${TCGA_BARCODE}_q1.bam
fi