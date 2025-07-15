#!/usr/bin/bash

STAR_bam_sorted_link=$1
TCGA_BARCODE=$2
TCGA_PROJECT=$3
THREADS=$4

# 1) Obtain the somatic SNVs VCF

sleep 10
VCF_snvs=$(readlink -f /g/strcombio/fsupek_cancer1/TCGA_bam/strelka_EVS/tissue_folders/${TCGA_PROJECT:5}/${TCGA_BARCODE}.vcf.gz)
slepp 10
STAR_bam_sorted_file=$(readlink $STAR_bam_sorted_link)

mkdir SNVs
cd SNVs

sleep 30
# If VCF is not found just abort mission buddy!
if [ -z "$VCF_snvs" ];
then 
	echo "somatic SNVS VCF not found"
else
        # 2) Normalize (left-align) VCF for indels
        conda activate bcftools
        gunzip -d -c $VCF_snvs > variants_somatic_snvs.vcf
        sleep 10
        bcftools norm -m - variants_somatic_snvs.vcf -f /local_scratch/gpalou/tmp/TCGA/GRCh38.d1.vd1.fa > variants_somatic_snvs_norm.vcf
        
        # 3) Tabix index VCF file
	conda activate /opt/anaconda3/envs/tabix
	bgzip -c variants_somatic_snvs_norm.vcf > strelka_somatic_snvs.vcf.gz
	tabix -p vcf strelka_somatic_snvs.vcf.gz       

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
                --forcedGT strelka_somatic_snvs.vcf.gz

        # 5) Run Strelka2 worflow
        ./StrelkaGermlineWorkflow/runWorkflow.py -m local -j ${THREADS}

        # 6) Copy output in current directory
        if [ -f "./StrelkaGermlineWorkflow/results/variants/variants.vcf.gz" ] 
        then
                mv ./StrelkaGermlineWorkflow/results/variants/variants.vcf.gz ./VCF_somatic_snvs_strelka_recall.vcf.gz
        fi
fi

# 2) Obtain the somatic INDELS VCF

sleep 10
VCF_indels=$(echo $VCF_snvs | sed 's/somatic.snvs.vcf.gz/somatic.indels.vcf.gz/')
sleep 10
mkdir indels
cd indels

sleep 30

if [ -z "$VCF_indels" ]
	then 
		echo "somatic INDELS VCF not found"
else
        # 2) Normalize (left-align) VCF for indels
        conda activate bcftools
        gunzip -d -c $VCF_indels > variants_somatic_indels.vcf
        sleep 10
        bcftools norm -m - variants_somatic_indels.vcf -f /local_scratch/gpalou/tmp/TCGA/GRCh38.d1.vd1.fa > variants_somatic_indels_norm.vcf
        
        # 3) Tabix index VCF file
	conda activate /opt/anaconda3/envs/tabix
	bgzip -c variants_somatic_indels_norm.vcf > strelka_somatic_indels.vcf.gz
	tabix -p vcf strelka_somatic_indels.vcf.gz

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
                --forcedGT strelka_somatic_indels.vcf.gz

        # 5) Run Strelka2 worflow
        ./StrelkaGermlineWorkflow/runWorkflow.py -m local -j ${THREADS}

        # 6) Copy output in current directory
        if [ -f "./StrelkaGermlineWorkflow/results/variants/variants.vcf.gz" ] 
        then
                mv ./StrelkaGermlineWorkflow/results/variants/variants.vcf.gz ../../VCF_somatic_indels_strelka_recall.vcf.gz
        fi
fi