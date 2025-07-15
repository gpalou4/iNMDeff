#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=ASE_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=20GB
##SBATCH --exclude=fsupeksvr5
#SBATCH --time=18:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/ASE/TCGA/"
tissue=$1 # TCGA-[X]
somatic_transcripts_filt=$2 # somatic_transcripts_filt // ""
VAF=$3 # e.g. 0.01
method=$4 # strelka // samtools
LOEUF_score=$5 # [0-9] // no

Rscript ${DIR}/ASE_NMD_efficiency.R ${DIR}/ASE_NMD_efficiency_PATHS_cluster.txt ${tissue} ${somatic_transcripts_filt} ${VAF} ${method} ${LOEUF_score}
