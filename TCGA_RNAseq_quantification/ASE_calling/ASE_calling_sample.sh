#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=ASE_calling
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --nodelist=fsupeksvr4
##SBATCH --mem-per-cpu=X
#SBATCH --mem=25GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
##SBATCH --exclude=fsupeksvr5
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

conda activate /home/dnaro/.conda/envs/NMD_regression_3

TCGA_CANCER=$1
TCGA_SAMPLE=$2
type=$3 # germline // somatic

PATH_FILE="/home/gpalou/projects/NMD/scripts/TCGA_RNAseq_quantification/ASE_calling"

Rscript ${PATH_FILE}/ASE_calling.R ${TCGA_CANCER} ${TCGA_SAMPLE} ${type}
