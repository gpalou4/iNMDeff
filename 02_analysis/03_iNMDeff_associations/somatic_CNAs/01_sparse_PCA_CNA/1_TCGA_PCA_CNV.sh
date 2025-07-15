#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=TCGA_pancancer_PCA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=65GB
##SBATCH --exclude=fsupeksvr2
##SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/PCA_CNV"

type=$1 # TCGA-[X] or pancancer
alpha=$2 # 0 to 0.0001
robust_SPCA=$3 # yes/no
num_PCs=$4

# Check somatic mut and CNV for that sample
#Rscript ${DIR}/1_TCGA_PCA_CNV.R ${type} ${alpha} ${robust_SPCA} ${num_PCs}
Rscript ${DIR}/1_TCGA_PCA_CNV.R ${type} ${alpha} ${robust_SPCA} ${num_PCs}

