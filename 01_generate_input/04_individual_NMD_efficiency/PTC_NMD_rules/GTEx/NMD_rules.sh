#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=NMD_rules_sample
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=50GB
#SBATCH --time=16:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/GTEx"
PATH_FILE=$1
GTEx_sample=$2
VCF_type=$3

conda env export

Rscript ${DIR}/NMD_rules.R ${DIR}/${PATH_FILE} ${GTEx_sample} ${VCF_type}



