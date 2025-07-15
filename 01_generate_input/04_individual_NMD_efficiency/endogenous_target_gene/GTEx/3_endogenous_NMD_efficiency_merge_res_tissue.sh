#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=merge_NMD_eff
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=6GB
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

WD="/home/gpalou/projects/NMD/scripts/NMD_efficiency/endogenous/ENSEMBL/GTEx/"
GTEx_tissue=$1
NMD_method=$2 #endogenous // ASE // PTC

Rscript ${WD}/3_endogenous_NMD_efficiency_merge_res_tissue.R ${GTEx_tissue} ${NMD_method}

