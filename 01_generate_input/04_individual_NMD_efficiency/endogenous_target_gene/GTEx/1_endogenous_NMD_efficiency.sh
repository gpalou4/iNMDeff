#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=GTEx_NMDeff_end
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=20GB
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

WD="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/GTEx"
GTEx_tissue=$1

Rscript ${WD}/1_endogenous_NMD_efficiency.R ${WD}/endogenous_NMD_efficiency_PATHS_cluster.txt ${GTEx_tissue} germline_transcripts_filt

