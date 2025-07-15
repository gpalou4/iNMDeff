#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=TCGA_samples_metadata
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=25GB
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

TCGA_cancer=$1
DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"

Rscript ${DIR}/TCGA_samples_metadata.R ${DIR}/NMD_rules_and_efficiency.txt ${TCGA_cancer}



