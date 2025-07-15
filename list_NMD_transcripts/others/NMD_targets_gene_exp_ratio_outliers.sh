#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=NMD_CL_outliers
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=60GB
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment

WD="/home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/TCGA"
PATHS=$1

Rscript ${WD}/NMD_targets_controls_outliers.R ${WD}/${PATHS}