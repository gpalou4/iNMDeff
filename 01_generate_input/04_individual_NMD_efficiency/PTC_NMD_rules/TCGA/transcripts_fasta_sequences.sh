#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=transcripts_fasta_sequences
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=12GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/gpalou/.conda/envs/gpalou_bedtools

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"

Rscript ${DIR}/transcripts_fasta_sequences.R



