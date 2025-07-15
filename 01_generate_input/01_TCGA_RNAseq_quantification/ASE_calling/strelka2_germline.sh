#!/usr/bin/bash

STAR_bam_sorted_link=$1
TCGA_BARCODE=$2
TCGA_PROJECT=$3
THREADS=$4

# 1) Obtain the germline VCF corresponding file
sleep 10
VCF=$(grep ${TCGA_BARCODE} /local_scratch/gpalou/tmp/TCGA_germline_variants/${TCGA_PROJECT:5}_VCFs_samples_list.txt)
sleep 10

STAR_bam_sorted_file=$(readlink $STAR_bam_sorted_link)

# conda activate samtools
# samtools sort $STAR_bam_sorted_file -o ${TCGA_BARCODE}.keep.merged.sorted.final.bam
# samtools index ${TCGA_BARCODE}.keep.merged.sorted.final.bam
# STAR_bam_sorted_file=${TCGA_BARCODE}.keep.merged.sorted.final.bam

sleep 30

# If VCF is not found just abort mission buddy!
if [ -z "$VCF" ];
then 
	echo "germline VCF not found"
else
        # 2) Normalize (left-align) VCF for indels
        conda activate bcftools
        gunzip -d -c $VCF > variants_germline.vcf
        sleep 10
        bcftools norm -m - variants_germline.vcf -f /local_scratch/gpalou/tmp/TCGA/GRCh38.d1.vd1.fa > variants_germline_norm.vcf

        # 3) Tabix index VCF file
	conda activate /opt/anaconda3/envs/tabix
	bgzip -c variants_germline_norm.vcf > strelka_germline.vcf.gz
	tabix -p vcf strelka_germline.vcf.gz

	# 3) Obtain the TCGA full barcode from the VCF header
	#VCF_header=$(zcat $VCF | grep "^#CHR")
	#TCGA_full_barcode=$(echo $VCF_header | cut -d ' ' -f10)

        # 4) Strelka2 configure file

        /opt/anaconda3/envs/strelka/bin/configureStrelkaGermlineWorkflow.py \
                --bam ${STAR_bam_sorted_file} \
                --referenceFasta /local_scratch/gpalou/tmp/TCGA/GRCh38.d1.vd1.fa \
                --exome \
                --rna \
                --callRegions /local_scratch/gpalou/tmp/chr.gz \
                --forcedGT strelka_germline.vcf.gz

        # 5) Run Strelka2 worflow
        ./StrelkaGermlineWorkflow/runWorkflow.py -m local -j ${THREADS}

        # 6) Copy output in current directory
        if [ -f "./StrelkaGermlineWorkflow/results/variants/variants.vcf.gz" ] 
        then
                mv ./StrelkaGermlineWorkflow/results/variants/variants.vcf.gz ./VCF_germline_strelka_recall.vcf.gz
        fi
fi