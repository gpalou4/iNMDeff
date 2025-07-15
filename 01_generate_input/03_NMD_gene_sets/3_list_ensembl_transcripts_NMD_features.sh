#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=list_ensembl_NMD_transcripts
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=150GB
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/gpalou/.conda/envs/gpalou_bedtools

WD="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/list_NMD_transcripts"
PATHS="ensembl_check_NMD_features_PATHS_cluster.txt"

Rscript ${WD}/3_list_ensembl_transcripts_NMD_features.R ${WD}/${PATHS}

