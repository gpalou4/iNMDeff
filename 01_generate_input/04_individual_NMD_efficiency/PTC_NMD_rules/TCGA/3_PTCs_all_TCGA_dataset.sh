#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=PTCs_all_TCGA_dataset
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --exclude=fsupeksvr5
#SBATCH --mem=12GB
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/PTC_NMD_rules/TCGA"
type=$1

Rscript ${DIR}/3_PTCs_all_TCGA_dataset.R ${type}



