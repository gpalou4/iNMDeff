#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=ASE_GTEx
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=16GB
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /g/strcombio/fsupek_home/dnaro/.conda/envs/NMD_regression_3

DIR="/g/strcombio/fsupek_home/gpalou/projects/NMD/scripts/NMD_efficiency/ASE/GTEx"
tissue=$1 # Liver
VAF=$2 # e.g. 0.01
LOEUF_score=$3 # [0-9] // no

Rscript ${DIR}/ASE_NMD_efficiency.R ${DIR}/ASE_NMD_efficiency_PATHS_cluster.txt ${tissue} ${VAF} ${LOEUF_score}


