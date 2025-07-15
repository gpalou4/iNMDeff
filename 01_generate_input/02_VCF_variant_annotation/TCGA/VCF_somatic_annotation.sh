#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=annovar_tissue_somatic
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=1GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

TCGA_project=$1

# Obtain list of samples from the tissue
TCGA_barcode_samples=$(ls /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA-${TCGA_project} | grep -v 'RNAseq')

# STAD
#TCGA_barcode_samples=$(ls /g/strcombio/fsupek_cancer2/TCGA_bam/STAD/ | grep "TCGA")
#TCGA_barcode_samples=$(ls /g/strcombio/fsupek_cancer2/TCGA_bam/OV/ | grep "TCGA" | grep -v "clinical")
#TCGA_barcode_samples=$(ls ~/data/TCGA_RNAseq_quantification/primary_tumor/TCGA-LAML/)

for TCGA_barcode in $TCGA_barcode_samples
do
        # SNVs
        VCF_somatic_snvs_sample=$(readlink -f /g/strcombio/fsupek_cancer1/TCGA_bam/strelka_EVS/tissue_folders/${TCGA_project}/${TCGA_barcode}.vcf.gz)
        if [ ! -f ${VCF_somatic_snvs_sample} ]; then
                echo "Sample doesn't have SNVs somatic VCF, skipping -->" ${VCF_somatic_snvs_sample}
                continue
        fi
        # Indels
        VCF_somatic_indels_sample=$(echo $VCF_somatic_snvs_sample | sed 's/somatic.snvs.vcf.gz/somatic.indels.vcf.gz/')
        # Annotation path
        VCF_somatic_annotation_path=$(echo $VCF_somatic_snvs_sample | sed "s/.*\/${TCGA_project}\(.*\)\/somatic.snvs.vcf.gz/\1/g")
        VCF_somatic_annotation_path="/g/strcombio/fsupek_cancer1/gpalou/TCGA_somatic_variants/"${TCGA_project}$VCF_somatic_annotation_path
        mkdir -p $VCF_somatic_annotation_path

	if [ -f "${VCF_somatic_annotation_path}/${TCGA_barcode}_annotated_snvs.HG_multianno.vcf.gz" ]; then
                echo "Already done, skipping -->" ${VCF_somatic_annotation_path}
                continue
        else
		sbatch VCF_somatic_annotation_sample.sh $TCGA_project $VCF_somatic_snvs_sample $VCF_somatic_indels_sample $VCF_somatic_annotation_path $TCGA_barcode
                #sleep 20
        fi
done

