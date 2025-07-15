#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=error_samples
##SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
##SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr
##SBATCH --mem-per-cpu=X
#SBATCH --mem=3GB
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

## It says how much time it will take the job to be scheduled
##SBATCH --test-only

TCGA_cancer=$1
DIR="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification"
GDC_TOKEN="gdc-user-token.2022-08-01T12_46_54.209Z.txt"
FASTQ_metadata="/g/strcombio/fsupek_home/gpalou/data/TCGA_fastq_RNAseq/TCGA_fastq_RNAseq_metadata.txt"

# Obtain samples
TCGA_samples=$(ls -lh /g/strcombio/fsupek_decider/gpalou/tmp/${TCGA_cancer} | grep -v ${TCGA_cancer} | awk '{print $9}')
#TCGA_samples=$(ls /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA-ACC/ | grep -v "TCGA-ACC")

for TCGA_sample in ${TCGA_samples}
    do
    # ASE Strelka2 output file for the sample
    ASE_output_file="/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_cancer}/${TCGA_sample}/ASE_germline/VCF_germline_strelka_recall.vcf.gz"
	# Obtain UUID
    UUID=$(grep $TCGA_sample $FASTQ_metadata | awk '{print $1}')
    # Move to FASTQ sample folder in DECIDER
    FASTQ_folder_path="/g/strcombio/fsupek_decider/gpalou/tmp/${TCGA_cancer}/${TCGA_sample}/fastq"
    cd ${FASTQ_folder_path}
    # Check if already re-downloaded
    month=$(ls -lh /g/strcombio/fsupek_decider/gpalou/tmp/${TCGA_cancer}/${TCGA_sample}/fastq | grep "fastq.tar.gz" | awk '{print $6}')
    # Download GDC FASTQ
    if [ -f ${ASE_output_file} ]; then
        echo "${TCGA_sample} sample already done, skipping -->" ${ASE_output_file}
    elif [ ${month} == "ago" ]; then
        echo "${TCGA_sample} sample already re-downloaded, skipping -->" ${FASTQ_folder_path}
    else
        echo "Downloading FASTQ files for --> ${TCGA_sample} at ${FASTQ_folder_path}"
        sbatch ${DIR}/ASE_calling/gdc_redownload.sh 4 ${DIR}/${GDC_TOKEN} 1 ${UUID}
        sleep 30
    fi
    done

