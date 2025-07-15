#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=endogenous_NB_regression
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=2GB
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

WD="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/TCGA"
PATHS_FILE=$1
TCGA_cancer=$2
TCGA_sample=$3
NMD_geneset=$4

Rscript ${WD}/endogenous_NMD_efficiency_regression.R ${WD}/${PATHS_FILE} ${TCGA_cancer} ${TCGA_sample} ${NMD_geneset}

