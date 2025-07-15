#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=NMD_features_ensembl
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --exclude=fsupeksvr5
#SBATCH --mem=40GB
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/gpalou/.conda/envs/gpalou_bedtools

WD="/home/gpalou/projects/NMD/scripts/list_NMD_transcripts"
PATHS=$1
N=$2
Rscript ${WD}/0_ensembl_check_NMD_features.R ${WD}/${PATHS} $N
