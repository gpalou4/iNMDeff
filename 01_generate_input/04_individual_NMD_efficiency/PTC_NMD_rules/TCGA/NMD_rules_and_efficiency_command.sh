#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=NMDrules_command
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr4
##SBATCH --mem-per-cpu=X
##SBATCH --mem=25GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
##SBATCH --exclude=fsupeksvr5
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

#conda activate /home/dnaro/.conda/envs/NMD_regression_3
PATH_FILE=$1
TCGA_cancer=$2
TCGA_sample=$3
VCF_type=$4
variant_type=$5
DIR="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"

if [ ${VCF_type} == "germline" ]; then
    TCGA_sample_output_file_1="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_cancer}/${TCGA_sample}/PTC_transcripts/AAgermline_PTC_transcripts.txt"
    TCGA_sample_output_file_2="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_cancer}/${TCGA_sample}/PTC_transcripts/AAgermline_PTC_transcripts_metadata.txt"
else
    TCGA_sample_output_file_1="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_cancer}/${TCGA_sample}/PTC_transcripts/AAsomatic_PTC_transcripts.txt"
    TCGA_sample_output_file_2="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_cancer}/${TCGA_sample}/PTC_transcripts/AAsomatic_PTC_transcripts_metadata.txt"
fi

if [ -f ${TCGA_sample_output_file_1} ]; then
        echo "${TCGA_sample} sample already done for PTC transcripts file, skipping -->" ${TCGA_sample_output_file_1}
    if [ -f ${TCGA_sample_output_file_2} ]; then
        echo "${TCGA_sample} sample already done for PTC transcripts metadata file, skipping -->" ${TCGA_sample_output_file_2}
    else
        # PTC metadata
        bash ${DIR}/PTCs_NMDeff_metadata.sh ${PATH_FILE} ${TCGA_cancer} ${TCGA_sample} ${VCF_type} ${variant_type}
    fi
else
    # Check NMD rules for all PTCs from the sample        
    bash ${DIR}/NMD_rules.sh ${PATH_FILE} ${TCGA_cancer} ${TCGA_sample} ${VCF_type} ${variant_type}
    bash ${DIR}/PTCs_NMDeff_metadata.sh ${PATH_FILE} ${TCGA_cancer} ${TCGA_sample} ${VCF_type} ${variant_type}
fi