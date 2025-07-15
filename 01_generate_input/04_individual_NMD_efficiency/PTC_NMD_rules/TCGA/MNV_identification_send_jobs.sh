#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=send_jobs
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=1GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

TCGA_cancer=$1
VCF_type=$2
variant_type=$3
DIR="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"

# Obtain samples
TCGA_samples=$(ls -lh /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/${TCGA_cancer}/ | grep -v ${TCGA_cancer} | awk '{print $9}')

for TCGA_sample in ${TCGA_samples}
    do
    echo ${TCGA_sample}
    sleep 1
    sbatch ${DIR}/MNV_identification.sh ${TCGA_cancer} ${TCGA_sample} ${VCF_type} ${variant_type}
    done




