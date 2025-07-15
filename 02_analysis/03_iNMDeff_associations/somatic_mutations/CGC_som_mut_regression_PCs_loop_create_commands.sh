#!/usr/bin/bash

#SBATCH --partition=normal_prio
##SBATCH --account=high_prio
#SBATCH --job-name=PCs_loop
#SBATCH --output=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.out
#SBATCH --error=/g/strcombio/fsupek_home/gpalou/out_err_files/%x-%j.err
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=X
##SBATCH --nodes=X
##SBATCH --mem-per-cpu=X
#SBATCH --mem=2GB
#SBATCH --time=3:00:00
#SBATCH --mail-type=ALL
##SBATCH --mail-user=guillermo.palou@irbbarcelona.org

# Activate base environment
conda activate /home/dnaro/.conda/envs/NMD_regression_3

DIR="/home/gpalou/projects/NMD/scripts/NMD_efficiency/analysis_results/somatic_associations"
OUT_DIR="/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/associations/somatic_associations/array_files"

type=$1 # Can be "cancers" or "pancancer"
NMD_method=$2 # Can be endogenous/PTC/ASE
test=$3 # Can be CGC/subtypes
VAF=$4 # VAF for the ASE method
alpha=$5 # can be one of the alphas tested
robust_SPCA=$6 # yes/no

if [[ ${NMD_method} == "ASE" ]]; then
    NMD_method_char=${NMD_method}_${VAF}
else
    NMD_method_char=${NMD_method}
fi

# Number of PCs for the cancer

> "${OUT_DIR}/${type}_${NMD_method_char}_${test}_som_mut_array_commands_${alpha}_robust_${robust_SPCA}.txt"

# Pancancer or cancer
if [[ ${type} == "pancancer" ]]; then
    num_PCs=1000
    # Loop through each PC
    for (( PC=1; PC<=$num_PCs; PC++ ))
    do
        echo "${DIR}/CGC_som_mut_regression_PCs_test.sh ${type} ${NMD_method} ${test} ${PC} ${VAF} ${alpha} ${robust_SPCA} " >> "${OUT_DIR}/${type}_${NMD_method_char}_${test}_som_mut_array_commands_${alpha}_robust_${robust_SPCA}.txt"
    done
else
    while read TCGA_cancer
    do  
        # This is, e.g., to use same PCA (but different sample PCs) from UCEC for both MSI and MSS TCGA_cancer
        TCGA_cancer_original=$(echo $TCGA_cancer | cut -d'_' -f1)
        num_PCs=$(head -1 /g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/${TCGA_cancer_original}_sparse_PCA_ind_${alpha}_robust_${robust_SPCA}.txt | wc -w)
        # Loop through each PC
        for (( PC=1; PC<=$num_PCs; PC++ ))
        do
            echo "${DIR}/CGC_som_mut_regression_PCs_test.sh ${TCGA_cancer} ${NMD_method} ${test} ${PC} ${VAF} ${alpha} ${robust_SPCA}" >> "${OUT_DIR}/${type}_${NMD_method_char}_${test}_som_mut_array_commands_${alpha}_robust_${robust_SPCA}.txt"
        done
    done < /g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt
fi








