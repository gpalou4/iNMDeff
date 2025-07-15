#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=PTCs_NMDeff_metadata
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
##SBATCH --mem=35GB
#SBATCH --time=16:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"
PATH_FILE=$1
TCGA_cancer=$2
TCGA_sample=$3
VCF_type=$4
variant_type=$5

Rscript ${DIR}/PTCs_NMDeff_metadata.R ${DIR}/${PATH_FILE} ${TCGA_cancer} ${TCGA_sample} ${VCF_type} ${variant_type}



