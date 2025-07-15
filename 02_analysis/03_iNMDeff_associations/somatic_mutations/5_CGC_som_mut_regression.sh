#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=CGC_regression_TCGA
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=20GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations"
TCGA_cancer=$1 # Can be TCGA-[X] or pancancer
NMD_method=$2 # Can be endogenous/PTC/ASE
test=$3 # Can be CGC/subtypes
CNV_PCs=$4 # 1 == Add PCs (the optimized ones), 0 == no
VAF=$5
alpha=$6 # can be one of the alphas tested
robust_SPCA=$7 # yes/no
PCA_PCs=$8 # num of PCs from sparse PCA (k)

# if [ ${NMD_method} == "ASE" ]; then
#     CNV_PCs=$(grep -w $TCGA_cancer /g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/${NMD_method}_${VAF}_best_PC_all_TCGA_cancers.txt | cut -f3)
# else
#     CNV_PCs=$(grep -w $TCGA_cancer /g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/optimizing_PCs/${NMD_method}_best_PC_all_TCGA_cancers.txt | cut -f3)
# fi

Rscript ${DIR}/5_CGC_som_mut_regression.R ${TCGA_cancer} ${NMD_method} ${test} ${CNV_PCs} ${VAF} ${alpha} ${robust_SPCA} ${PCA_PCs}
