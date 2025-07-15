#!/usr/bin/bash

#SBATCH --partition=normal_prio
#SBATCH --job-name=cNMDeff
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
##SBATCH --nodelist=fsupeksvr2
##SBATCH --excludenode=fsupeksvr5
#SBATCH --mem=4GB
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

NMD_gene=$1

Rscript /home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/Estimate_cNMDeff/Estimate_cNMDeff.R ${NMD_gene}

