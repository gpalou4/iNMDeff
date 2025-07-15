#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=NMD_efficiency_endogenous
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=12GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

WD="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/TCGA"
PATHS_FILE=$1
TCGA_tissue=$2

Rscript ${WD}/endogenous_NMD_efficiency.R ${WD}/${PATHS_FILE} ${TCGA_tissue} somatic_transcripts_filt germline_transcripts_filt

