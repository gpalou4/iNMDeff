#!/usr/bin/bash

#SBATCH --partition=normal_prio_long
##SBATCH --account=high_prio
#SBATCH --job-name=CGC_regression_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=12GB
##SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations"
TCGA_cancer=$1 # Can be TCGA-[X] or pancancer
NMD_method=$2 # Can be endogenous/PTC/ASE
test=$3 # Can be CGC/subtypes
CNV_PCs=$4 # Can be any number
VAF=$5 # VAF for the method
alpha=$6 # can be one of the alphas tested
robust_SPCA=$7 # yes/no

Rscript ${DIR}/CGC_som_mut_regression_PCs_test.R ${TCGA_cancer} ${NMD_method} ${test} ${CNV_PCs} ${VAF} ${alpha} ${robust_SPCA}

